# Install relevant library for HTTP requests
library(httr)

# Set gene_id variable
variantId <- "17_44352876_C_T"

# Build query string
query_string = "
query v2g($variantId: String!) {
  genesForVariant(variantId: $variantId) {
    gene {
      id
    }
    variant
    overallScore
    distances {
      sourceId
      aggregatedScore
      tissues {
      	distance
      }
    }
  }
}"

# Set base URL of GraphQL API endpoint
base_url <- "https://api.genetics.opentargets.org/graphql"

# Set variables object of arguments to be passed to endpoint
variables <- list("variantId" = variantId)

# Construct POST request body object with query string and variables
post_body <- list(query = query_string, variables = variables)

# Perform POST request
r <- POST(url=base_url, body=post_body, encode='json')

df = content(r)
# Print first entry of V2G data console
head(content(r)$data$genesForVariant, 1)

# Flatten the nested result fields into a dataframe
library(rlist)
list_result = content(r)$data$genesForVariant
x = lapply(list_result, list.flatten)

library(dplyr)
df = bind_rows(x)