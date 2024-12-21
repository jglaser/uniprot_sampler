from SPARQLWrapper import SPARQLWrapper, JSON
import pandas as pd
import sys
import concurrent.futures
import queue
from retrying import retry

endpoint = "https://sparql.uniprot.org/sparql"
sparql = SPARQLWrapper(endpoint)
sparql.setReturnFormat(JSON)
sparql.setTimeout(60)

NUM_WORKERS=16
identity = sys.argv[1]

query_template = """
PREFIX up: <http://purl.uniprot.org/core/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX faldo: <http://biohackathon.org/resource/faldo#>

SELECT DISTINCT
    (SUBSTR(STR(?protein), 33) AS ?uniprot)
    (GROUP_CONCAT(DISTINCT STR(?ec_number); separator=", ") AS ?ec_numbers)
    (GROUP_CONCAT(DISTINCT CONCAT(?go_number, " (", STR(?go_term), ")"); separator=", ") AS ?go_terms)
    (GROUP_CONCAT(DISTINCT CONCAT(
        "Position: ", COALESCE(STR(?binding_start), "-"), "-", COALESCE(STR(?binding_end), "-"),
        "; Ligand: ", COALESCE(STR(?chebi_label), "N/A")
    ); separator=", ") AS ?binding_sites)
    (GROUP_CONCAT(DISTINCT CONCAT(
        "Position: ", COALESCE(STR(?active_start), "-"), "-", COALESCE(STR(?active_end), "-"),
        "; Description: ", COALESCE(STR(?active_desc), "N/A")
    ); separator=", ") AS ?active_sites)
    (GROUP_CONCAT(DISTINCT CONCAT(
        "Position: ", COALESCE(STR(?domain_start), "-"), "-", COALESCE(STR(?domain_end), "-"),
        "; Description: ", COALESCE(STR(?domain_desc), "N/A")
    ); separator=", ") AS ?domains)
    (GROUP_CONCAT(DISTINCT CONCAT(
        "Cluster: ", STR(?cluster), "; Identity: ", STR(?identity)
    ); separator=", ") AS ?uniref_clusters)
    ?protein_name
    (GROUP_CONCAT(DISTINCT STR(?function_desc); separator=", ") AS ?protein_function)
    ?sequence
WHERE {
  {
    SELECT DISTINCT ?cluster ?protein ?ec_number
    WHERE {
        {
            # Inner query: Filter clusters with identity = 0.5 and retrieve proteins
            SELECT DISTINCT ?cluster ?protein
            WHERE {
                ?cluster up:identity ?identity .
                FILTER(?identity = _identity_)
                FILTER ( 1 >  <SHORT_OR_LONG::bif:rnd> (_decimate_, ?cluster, ?protein))

                ?protein up:representativeFor ?cluster .
            }
        }

        # Outer query: Retrieve EC numbers for the proteins
        OPTIONAL {
            ?protein (up:enzyme | up:domain/up:enzyme | up:component/up:enzyme) ?enzyme .
            FILTER(STRSTARTS(STR(?enzyme), "http://purl.uniprot.org/enzyme/"))
            BIND(STRAFTER(STR(?enzyme), "http://purl.uniprot.org/enzyme/") AS ?ec_number)
        }

        # Ensure proteins with EC numbers are included
        FILTER(?ec_number != "")
    }
  }
  # Retrieve sequence
  OPTIONAL {
    ?protein up:sequence ?sequence_obj .
    ?sequence_obj rdf:value ?sequence .
  }

  # Retrieve protein name
  OPTIONAL {
    ?protein rdfs:label ?protein_name .
  }

  # Retrieve active sites
  OPTIONAL {
    ?protein up:annotation ?active_annotation .
    ?active_annotation a up:Active_Site_Annotation .
    ?active_annotation up:range ?active_range .
    OPTIONAL {
        ?active_range faldo:begin/faldo:position ?active_start .
        ?active_range faldo:end/faldo:position ?active_end .
    }
    OPTIONAL {
        ?active_annotation rdfs:comment ?active_desc .
    }
  }
    # Retrieve GO terms with numbers
    OPTIONAL {
        ?protein up:classifiedWith ?go_class .
        ?go_class rdfs:label ?go_term .
        BIND(STRAFTER(STR(?go_class), "http://purl.obolibrary.org/obo/GO_") AS ?go_number)
    }

    # Retrieve binding sites with ChEBI ligand names
    OPTIONAL {
        ?protein up:annotation ?binding_annotation .
        ?binding_annotation a up:Binding_Site_Annotation .
        ?binding_annotation up:range ?binding_range .
        OPTIONAL {
            ?binding_range faldo:begin/faldo:position ?binding_start .
            ?binding_range faldo:end/faldo:position ?binding_end .
        }
        OPTIONAL {
            ?binding_annotation up:ligand ?ligand .
            ?ligand rdfs:subClassOf ?chebi_uri .
            OPTIONAL {
                ?chebi_uri rdfs:label ?chebi_label .
            }
        }
    }
    # Retrieve domains (from Region_Annotation)
    OPTIONAL {
        ?protein up:annotation ?region_annotation .
        ?region_annotation a up:Region_Annotation .
        OPTIONAL {
            ?region_annotation up:range ?region_range .
            ?region_range faldo:begin/faldo:position ?domain_start .
            ?region_range faldo:end/faldo:position ?domain_end .
        }
        OPTIONAL {
            ?region_annotation rdfs:comment ?domain_desc .
        }
    }

    # Retrieve function annotation
    OPTIONAL {
        ?protein up:annotation ?function_annotation .
        ?function_annotation a up:Function_Annotation .
        OPTIONAL {
            ?function_annotation rdfs:comment ?function_desc .
        }
    }

}
GROUP BY ?protein ?sequence ?protein_name
ORDER BY ?uniprot
LIMIT 1000
"""

@retry(stop_max_attempt_number=3,
       stop_max_delay=5000, # 5s
       )
def fetch(decimate):
    query = query_template.replace('_identity_',identity)
    query = query.replace('_decimate_',f'{decimate}')

    sparql.setQuery(query)
    results = sparql.query().convert()
    entries = []
    for result in results["results"]["bindings"]:
        entry = {
            "UniProt ID": result["uniprot"]["value"],
            "EC Numbers": result["ec_numbers"]["value"] if "ec_numbers" in result else "",
            "GO Terms": result["go_terms"]["value"] if "go_terms" in result else "",
            "Binding Sites": result["binding_sites"]["value"] if "binding_sites" in result else "",
            "Active Sites": result["active_sites"]["value"] if "active_sites" in result else "",
            "Domains": result["domains"]["value"] if "domains" in result else "",
            "UniRef Clusters": result["uniref_clusters"]["value"] if "uniref_clusters" in result else "",
            "Protein Name": result["protein_name"]["value"] if "protein_name" in result else "",
            "Function": result["protein_function"]["value"] if "protein_function" in result else "",
            "Sequence": result["sequence"]["value"] if "sequence" in result else "",
        }
        entries.append(entry)
    return pd.DataFrame(entries)

# Execute
if __name__ == "__main__":
    n_batches = int(sys.argv[2])
    decimate = int(sys.argv[3])

    from tqdm.contrib.concurrent import thread_map

    # fetch n random batches
    results = thread_map(lambda x: fetch(decimate),
                         range(n_batches),
                         max_workers=NUM_WORKERS)
    df = pd.concat(results)
    df = df.drop_duplicates(subset='UniProt ID').reset_index(drop=True)
    fn = f'uniref_{identity}_{n_batches}_{decimate}.csv'
    df.to_csv(fn)
    print(f'Fetched {len(df)} unique records and save as {fn}.')
