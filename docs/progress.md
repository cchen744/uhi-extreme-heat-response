2026.02.27
What I am currently doing is to get FeatureCollection where each feature is one cell row with properties I can export.
The ideal data structure looks like:
| cell_id| suhi | delta_uhi | (optional built vars: lcz, impervious, ndvi, …) |

So the idea is, instead of looking at UHI at the city level, I want to break the city into smaller cells. And for each cell, I extract daily temperature and compute  delta_uhi.
This delta-uhi indicates how extreme heat intensifies uhi compared to normal days.After that I will connect the table with LCZ to further analyze the
influence of built-up environment.

My current technical issue is:
1. computational power seems to be wasted too much - even a test sample could spend over half an hour.
2. Since I have to go through several filtering process - there are multiple sequential steps: spatial masking, QC filtering, temporal aggregation, and then joining across datasets.
This one really trips me over because  I find even one tiny mismatch in column names, a null value, different scale in pandas an gee...could lead to an empty output.
So a lot of my time right now is spent on debugging data consistency across these steps rather than the analysis itself. Hopefully next week I can start on the LCZ join.

这个部分代码逻辑有问题：
                urb.reduceRegions(
                    collection=grid,
                    reducer=combined,
                    scale=lst_scale_m,
                    crs=crs,
                    tileScale=tileScale
                )
                .filter(ee.Filter.notNull([key_val]))
                .map(lambda ft: ft.set({
                    "date":      date_str,
                    "SUHI_cell": ee.Number(ft.get(key_val)).subtract(ee.Number(lst_rur)),
                    "cell_n":    ft.get(key_cnt),
                    "rural_n":   rural_n,
                }).select(["date", "cell_id", "SUHI_cell", "cell_n", "rural_n"]))
            )

            # 没有识别出extreme_day然后做aggregated_extremheat_uhi - aggregated_normla_uhi
