import pandas as pd
import requests
import json
import databases


url = "http://l2s2.maayanlab.cloud/graphql"

def enrich_l2s2_single_set(geneset: list, db, first=1000):
    query = {
    "operationName": "EnrichmentQuery",
    "variables": {
        "filterTerm": " ",
        "offset": 0,
        "first": first,
        "filterFda": False,
        "sortBy": "pvalue_up",
        "filterKo": False,
        "genes": geneset,
    },
    "query": """query EnrichmentQuery(
                    $genes: [String]!
                    $filterTerm: String = ""
                    $offset: Int = 0
                    $first: Int = 10
                    $filterFda: Boolean = false
                    $sortBy: String = ""
                    $filterKo: Boolean = false
                    ) {
                    currentBackground {
                        enrich(
                        genes: $genes
                        filterTerm: $filterTerm
                        offset: $offset
                        first: $first
                        filterFda: $filterFda
                        sortby: $sortBy
                        filterKo: $filterKo
                        ) {
                        nodes {
                            geneSetHash
                            pvalue
                            adjPvalue
                            oddsRatio
                            nOverlap
                            geneSets {
                            nodes {
                                term
                                id
                                nGeneIds
                                geneSetFdaCountsById {
                                nodes {
                                    approved
                                    count
                                }
                                }
                            }
                            totalCount
                            }
                        }
                        totalCount
                        consensusCount
                        consensus {
                            drug
                            oddsRatio
                            pvalue
                            adjPvalue
                            approved
                            countSignificant
                            countInsignificant
                            countUpSignificant
                            pvalueUp
                            adjPvalueUp
                            oddsRatioUp
                            pvalueDown
                            adjPvalueDown
                            oddsRatioDown
                        }
                        }
                    }
                    }
                    """,
    }

    headers = {
        "Accept": "application/json",
        "Content-Type": "application/json"
    }

    response = requests.post(url, data=json.dumps(query), headers=headers)

    response.raise_for_status()
    res = response.json()

    #consensus = pd.DataFrame(res['data']['currentBackground']['enrich']['consensus'])
    consensus = res['data']['currentBackground']['enrich']['consensus']
    #enrichment = pd.DataFrame(res['data']['currentBackground']['enrich']['nodes'])
    enrichment = res['data']['currentBackground']['enrich']['nodes']# %%
    df_consensus = pd.DataFrame(consensus).rename(columns={'drug': 'perturbation'})

    df_enrichment = pd.json_normalize(
        enrichment, 
        record_path=['geneSets', 'nodes'], 
        meta=['geneSetHash', 'pvalue', 'adjPvalue', 'oddsRatio', 'nOverlap']
    )
    if df_enrichment.empty:
        return pd.DataFrame(), pd.DataFrame()
    # df_enrichment["approved"] = df_enrichment["geneSetFdaCountsById.nodes"].map(lambda x: x[0]['approved'] if len(x) > 0 else False)
    # df_enrichment["count"] = df_enrichment["geneSetFdaCountsById.nodes"].map(lambda x: x[0]['count'] if len(x) > 0 else 0)
    # df_enrichment.drop(columns=['geneSetFdaCountsById.nodes'], inplace=True)
    # df_enrichment['batch'] = df_enrichment["term"].map(lambda t: t.split('_')[0])
    # df_enrichment["timepoint"] = df_enrichment["term"].map(lambda t: t.split('_')[1])
    df_enrichment["Cell_Line"] = df_enrichment["term"].map(lambda t: t.split('_')[2])
    # df_enrichment["batch2"] = df_enrichment["term"].map(lambda t: t.split('_')[3])
    
    df_enrichment["DrugBank_Name"] = df_enrichment["term"].map(lambda t: pd.NA if len(t.split('_')[4].split(' ')) == 2 else t.split('_')[4])
    
    # df_enrichment['concentration'] = df_enrichment["term"].map(lambda t: t.split('_')[5].split(' ')[0] if len(t.split('_')) > 5 else "N/A")
    # df_enrichment['direction'] = df_enrichment["term"].map(lambda t: t.split(' ')[1])\

    df_enrichment = df_enrichment.dropna()

    def name2DBid(name):
        try:
            return db.search(name).id
        except:
            # print(name)
            return pd.NA

    df_enrichment["DrugBank_ID"] = df_enrichment["DrugBank_Name"].apply(name2DBid)

    df_enrichment = df_enrichment.dropna()

    return df_enrichment

db = databases.DrugBank()

data = enrich_l2s2_single_set([
      "KIF1B",
      "ST3GAL3",
      "PLPP3",
      "JAK1",
      "PDE4B",
      "AHCYL1",
      "PDE4DIP",
      "DISP1",
      "AC119674.1"
    ], db)

print(data)