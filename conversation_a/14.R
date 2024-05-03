# Compute proportion of sgCntrl cells within each cluster
cluster_summary <- tsne_data %>%
  group_by(cluster, group_id) %>%
  summarise(count = n()) %>%
  spread(group_id, count, fill = 0) %>%
  mutate(total = rowSums(select(., -cluster)),
         sgCntrl_prop = sgCntrl / total) %>%
  arrange(desc(sgCntrl_prop))

# View the cluster summary
print(cluster_summary)

