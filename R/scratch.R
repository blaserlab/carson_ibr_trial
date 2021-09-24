

colData(cds_aligned) %>%
  as_tibble() %>%
  filter(partition_assignment_1 %in% c("T", "pyr-T")) %>%
  # filter(treatment == "ibrutinib") %>%
  group_by(patient, response, treatment, partition_assignment_1) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = partition_assignment_1, values_from = n) %>%
  mutate(pyr_ratio = `pyr-T`/`T`) %>%
  mutate(log2pyr_ratio = log2(pyr_ratio)) %>%
  ggplot(mapping = aes(x = response, y = log2pyr_ratio, fill = response)) +
  geom_jitter(
    shape = jitter_shape,
    width = jitter_width,
    size = jitter_size,
    stroke = jitter_stroke,
    height = 0
  ) +
  stat_summary(
    fun.data = data_summary_mean_se,
    color = summarybox_color,
    size = summarybox_size,
    width = summarybox_width,
    alpha = summarybox_alpha,
    geom = summarybox_geom
  ) +
  theme(legend.position = "none")+
  facet_wrap(facets = vars(treatment)) +
  theme(strip.background = element_blank())
  


glm_data <- colData(cds_aligned) %>%
  as_tibble() %>%
  filter(partition_assignment_1 %in% c("T", "pyr-T")) %>%
  # filter(treatment == "ibrutinib") %>%
  group_by(patient, response, disease, treatment, partition_assignment_1) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = partition_assignment_1, values_from = n) %>%
  mutate(pyr_ratio = `pyr-T`/`T`) %>%
  mutate(log2pyr_ratio = log2(pyr_ratio))
glm_data


glm_res1 <- glm(log2pyr_ratio ~ response, data = glm_data, family = "gaussian")
summary(glm_res1)

glm_res2 <- glm(log2pyr_ratio ~ response+disease, data = glm_data, family = "gaussian")
summary(glm_res2)

glm_res3 <- glm(log2pyr_ratio ~ response+treatment, data = glm_data, family = "gaussian")
summary(glm_res3)


