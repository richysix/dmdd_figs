#!/usr/bin/env python

# This is an automatically generated script to run your query
# to use it you will require the intermine python client.
# To install the client, run the following command from a terminal:
#
#     sudo easy_install intermine
#
# For further documentation you can visit:
#     http://intermine.readthedocs.org/en/latest/web-services/

# This script has been edited from the original automatically generated one
# to alter the constraints and columns returned.
# It was downloaded from http://www.mousemine.org/mousemine/template.do?name=Gene_Expression&scope=all

# The following two lines will be needed in every python script:
from intermine.webservice import Service
service = Service("http://www.mousemine.org/mousemine/service")

# Get a new query on the class (table) you will be querying:
query = service.new_query("GXDExpression")

# The view specifies the output columns
query.add_view(
    "assayType", "feature.symbol", "feature.primaryIdentifier", "stage", "age",
    "structure.name", "structure.identifier", "emaps", "strength", "pattern",
    "genotype.symbol", "sex", "assayId", "probe",
    "image", "publication.mgiJnum", "publication.pubMedId"
)

# This query's custom sort order is specified below:
query.add_sort_order("GXDExpression.assayId", "ASC")

# You can edit the constraint values below
query.add_constraint("feature.organism.taxonId", "=", "10090", code = "B")
#query.add_constraint("feature", "LOOKUP", "Pax6", code = "A")

# Uncomment and edit the code below to specify your own custom logic:
# query.set_logic("A and B")

# header
print "\t".join( [ "assayType", "gene.name", "MGI.id", \
    "stage", "age", "structure.name", "EMAPA", "EMAPS", "strength", "pattern", \
    "genotype.symbol", "sex", "mgiJnum", "pubMedId" ] )

for row in query.rows():
    print "\t".join( str(x) for x in [ row["assayType"], row["feature.symbol"], \
        row["feature.primaryIdentifier"], row["stage"], row["age"],
        row["structure.name"], row["structure.identifier"], row["emaps"], \
        row["strength"], row["pattern"], \
        row["genotype.symbol"], row["sex"], \
        row["publication.mgiJnum"], row["publication.pubMedId"] ] )

