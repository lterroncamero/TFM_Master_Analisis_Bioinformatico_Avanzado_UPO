
library(rbioapi)
library(STRINGdb)
library(dplyr)

# Conexion
string_db <- STRINGdb$new( version="11.5", species=9606, score_threshold=500, input_directory="")


# load_data
data_endo=read.delim("Emdometrio.tsv", dec = ",")[,c(2,3,4,6,7)]
data_endo=data_endo %>% filter(adj.P.Val < 0.05)
data_embryo=read.delim("1.Embrion.tsv", dec = ",")[,c(-1,-3,-6)]
data_embryo=data_embryo %>% filter(FDR < 0.05)
colnames(data_endo)=colnames(data_embryo)

data_endo$tissue= "Endometrium"
data_embryo$tissue= "Embryo"

data=rbind(data_endo, data_embryo)

# map
data_mapped = string_db$map( data, "Gene", removeUnmappedRows = TRUE )
head(data_mapped)

# Write file for cytoscape
duplicados = data_mapped[duplicated(data_mapped$STRING_id),"STRING_id"]
data_mapped[data_mapped$STRING_id %in% duplicados,"tissue"]="Endometrium/Embryo"

data_mapped = data_mapped %>% distinct(STRING_id, .keep_all = T) %>% relocate(STRING_id, .before = Gene)

write.table(data_mapped, "File.tsv", sep="\t", row.names = F, quote = F)


# get interactions

interactions= rba_string_interactions_network(data_mapped$STRING_id, species = 9606, required_score = 500)
head(interactions)

# nscore: gene neighborhood score
# fscore: gene fusion score
# pscore: phylogenetic_profile_score
# ascore: "co-expression_score"
# escore: "experimental_score"
# dscore: "database_score"
# tscore: "textmining_score"


# remove textmining
interactions_filt = unique(interactions) %>% filter((nscore+fscore+pscore+ascore+escore+dscore) > 0)
# Actualizar score
interactions_filt$score = rowSums(interactions_filt[,c("nscore","fscore","pscore","ascore","escore","dscore")])
# Añadir num de evidencias de cada interaccion
interactions_filt$n_suport= apply(interactions_filt[,c("nscore","fscore","pscore","ascore","escore","dscore")], 1, FUN = function(x) sum(x>0))
interactions_filt=interactions_filt[,-c(3,4,5,12)]
dim(interactions_filt)

colnames(interactions_filt)<-c("stringdb_A", "stringdb_B", "score", "gene_neighborhood_score", "gene_fusion_score", "phylogenetic_profile_score", "co-expression_score", "experimental_score", "textmining_score","n_suport")
write.table(interactions_filt, "Edges.tsv", sep="\t", row.names = F, quote = F)

# Color
data_mapped$color = "#FC0000"
data_mapped[data_mapped$tissue=="Embryo", "color"]="#0400FC"
# post payload information to the STRING server
payload_id <- string_db$post_payload( data_mapped$STRING_id, colors=data_mapped$color )

# display a STRING network png with the "halo"
string_db$plot_network( data_mapped$STRING_id, payload_id=payload_id )

# grafo and link
data_graph= string_db$get_subnetwork(data_mapped$STRING_id)  # returns a subgraph from the given input proteins

string_db$get_link(data_mapped$STRING_id, required_score=500, network_flavor="evidence", payload_id = payload_id)
data_inter= string_db$get_interactions(data_mapped$STRING_id)
