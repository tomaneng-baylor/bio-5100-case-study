# Compute proportion of sgCntrl cells within each cluster
cluster_summary <- tsne_data %>%
  group_by(cluster, group_id) %>%
  summarise(count = n()) %>%
  spread(group_id, count, fill = 0) %>%
  ungroup() %>%
  mutate(total = rowSums(select(., -cluster)),
         sgCntrl_prop = ifelse(total == 0, 0, `sgCntrl` / total)) %>%
  arrange(desc(sgCntrl_prop))

