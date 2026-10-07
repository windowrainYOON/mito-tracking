# Mito Analyzer 지표 설명

결과 표(`_per_cell.csv`, `_per_mito.csv`, `_results.xlsx`)와 *Analysis* 탭의 히트맵·그래프에 나오는 지표를 설명합니다.
모든 길이는 µm, 넓이는 µm² 입니다. green 세기는 배경(세포 사이 어두운 틈의 밝기)을 뺀 값입니다.

용어
- **미토 마스크 (footprint)**: 세포 안에서 미토콘드리아로 판정된 픽셀 전체입니다. 기본 설정(*Split objects*)에서는 국소 임계값으로 잡고, 붙어 있는 미토콘드리아를 어둡고 좁은 접촉선에서 나눈 결과입니다.
- **골격 (skeleton)**: 미토 마스크를 1픽셀 두께 선으로 줄인 것입니다. MiNA(Fiji)와 같은 방식으로 가지(branch), 분기점(junction), 끝점(end)을 셉니다.
- **네트워크 (network)**: 골격에서 서로 연결된 덩어리 하나입니다.
- **미토 조각 (mito object)**: 미토 마스크에서 서로 떨어진 조각 하나입니다 (0.05 µm² 미만은 제외).
- **green puncta**: green 이미지에서 주변보다 뚜렷하게 밝은 점(top-hat 후 이미지마다 하나의 임계값).

---

## 1. 미토콘드리아 형태 (히트맵의 가로축, X)

### 양 (얼마나 많은가)
| 표시 이름 | 열 이름 | 뜻 | 해석 |
|---|---|---|---|
| Mito footprint fraction | `footprint_fraction` | 미토 마스크 넓이 ÷ 세포 넓이 | 세포가 미토콘드리아로 얼마나 채워져 있는지. 세포 크기와 무관한 양의 지표 |
| Mito footprint (µm²) | `mitochondrial_footprint_um2` | 미토 마스크 넓이 | 미토콘드리아 총량(면적). 큰 세포일수록 커짐 |
| Total branch length (µm) | `total_branch_length_um` | 골격 가지 길이의 합 | 미토콘드리아 총 길이 |

### 연결성 (네트워크가 얼마나 이어져 있는가) — MiNA 지표
| 표시 이름 | 열 이름 | 뜻 | 해석 |
|---|---|---|---|
| Networks | `n_networks` | 골격 덩어리 수 | 많을수록 잘게 나뉨 (세포 크기에도 비례) |
| Branches | `n_branches` | 골격 가지 수 | 미토 총량과 분기를 함께 반영 |
| Branch length mean (µm) | `branch_length_mean_um` | 가지 하나의 평균 길이 | 길수록 길쭉한 튜브형, 짧을수록 짧은 조각·점 |
| Network length mean (µm) | `summed_branch_lengths_mean_um` | 네트워크 하나에 속한 가지 길이 합의 평균 | 네트워크 하나의 평균 크기. 융합(fusion)이 많으면 커짐 |
| Branches per network | `network_branches_mean` | 네트워크당 평균 가지 수 | 1에 가까우면 가지 없는 막대·점, 클수록 그물 모양 |
| Donuts | `donuts` | 고리 모양으로 닫힌 네트워크 수 | 고리형(링) 미토콘드리아. 스트레스·부종 때 늘기도 함 |

### 조각화 (얼마나 잘게 쪼개졌는가) — 미토 조각 단위
| 표시 이름 | 열 이름 | 뜻 | 해석 |
|---|---|---|---|
| Mito objects | `n_mito_objects` | 세포 안 미토 조각 수 | 세포 크기에 비례하므로 아래 밀도 지표와 함께 볼 것 |
| Mito object length mean (µm) | `mito_length_mean_um` | 조각 하나의 골격 길이 평균 | 길수록 융합·튜브형, 짧을수록 분열(fission) |
| Mito aspect ratio mean | `mito_aspect_ratio_mean` | 조각의 장축 ÷ 단축 평균 | 1 = 둥근 점, 클수록 길쭉 |
| Objects per 100 µm² footprint | `objects_per_100um2_footprint` | 미토 면적 100 µm² 당 조각 수 | **대표적인 조각화 지표.** 클수록 잘게 쪼개짐 |
| Area-weighted object size (µm²) | `area_weighted_mean_object_um2` | Σ(조각 넓이²) ÷ Σ(조각 넓이) | "미토 픽셀 하나가 평균적으로 속한 조각의 크기". 작을수록 조각화. 작은 부스러기 수에 덜 흔들림 |
| Form factor mean | `form_factor_mean` | 조각별 둘레² ÷ (4π × 넓이)의 평균 | 1 = 원, 클수록 길고 구불구불하거나 가지가 많음 |
| Footprint in small round objects | `frac_footprint_small_round` | 미토 면적 중 작고(< 0.5 µm²) 둥근(장축/단축 < 2) 조각이 차지하는 비율 | **분열된 점 모양 미토의 비율.** 클수록 조각화 |

조각화가 심해지면 보통 함께 움직이는 방향:
`Objects per 100 µm²` ↑, `Footprint in small round objects` ↑ / `Branch length mean`, `Network length mean`, `Branches per network`, `Mito object length mean`, `Form factor mean`, `Area-weighted object size` ↓.
그래서 히트맵에서 이 두 묶음의 색이 서로 반대로 나오는 것이 정상입니다.

---

## 2. Green(관심 단백질) 신호 (히트맵의 세로축, Y)

### 세기와 분포
| 표시 이름 | 열 이름 | 뜻 | 해석 |
|---|---|---|---|
| Green mean, whole cell | `green_poi_mean_cell` | 세포 전체 green 평균 | 발현량 |
| Green mean on mito | `green_mean_on_mito` | 미토 마스크 위 green 평균 | 미토콘드리아에 있는 단백질 농도 |
| Green mean off mito | `green_mean_off_mito` | 미토 밖(세포질·핵) green 평균 | 미토에 가지 못하고 남은 신호 |
| Green enrichment on/off mito | `green_enrichment_on_mito` | on ÷ off | 미토로 얼마나 농축되었나 (발현량과 무관한 비율) |
| Fraction of green on mito | `green_fraction_on_mito` | 미토 위 green 합 ÷ 세포 green 합 | 세포 green의 몇 %가 미토에 있나. **표적화(targeting) 효율** |
| Green integrated on mito | `green_integrated_on_mito` | 미토 위 green 합 × 픽셀 넓이 | 미토에 있는 단백질 총량 (미토 양과 발현량 모두 반영) |
| Pearson red–green | `pearson_red_green` | 세포 안 모든 픽셀에서 red와 green의 Pearson 상관 | 1에 가까울수록 green이 미토 모양을 그대로 따라감 (고른 분포) |
| Bright green area (% cytoplasm) | `green_bright_area_percent` | 세포질(세포 − 핵) 중 밝은 green이 덮은 비율 | green+/− 판정에 쓰는 값 (2 % 이상 = green+) |

### Puncta (밝은 점) 기반
| 표시 이름 | 열 이름 | 뜻 | 해석 |
|---|---|---|---|
| Green puncta | `n_green_puncta` | 세포 안 green 점 수 | 세포 크기에 비례 |
| Green puncta on mito | `n_puncta_on_mito` | 미토와 겹치는 점 수 | |
| Fraction of puncta on mito | `fraction_puncta_on_mito` | 미토와 겹치는 점 ÷ 전체 점 | 점이 미토 위에 생기나, 밖에 생기나 |
| Puncta per 100 µm² | `green_puncta_density_per_100um2` | 세포 넓이 100 µm² 당 점 수 | 점(응집·foci)이 얼마나 많나 |
| Manders: puncta green on mito | `manders_m_green` | 점 안 green 중 미토 위에 있는 비율 | 점의 green이 미토에 있는 정도 |
| Manders: mito red under puncta | `manders_m_mito` | 미토 red 중 점과 겹치는 비율 | **미토 중 얼마만큼이 green 점으로 덮였나.** green이 미토 전체에 고르게 퍼지면 크고, 일부 점에만 모이면 작음 |
| Fraction of mito objects with puncta | `fraction_mito_with_puncta` | green 점이 하나라도 겹치는 미토 조각의 비율 | 미토 조각 중 몇 %가 green을 가졌나 |

### 미토 조각 단위 지표 (*Per mito object* 레벨)
X: `length_um`(골격 길이), `aspect_ratio`, `area_um2`, `n_branches`, `n_junctions`, `n_endpoints`, `solidity`(볼록 껍질 대비 채움 정도, 1 = 꽉 참), `form_factor`, `major_axis_um`, `red_mean`.
Y: `green_mean`, `green_integrated`, `green_max`, `n_green_puncta`, `puncta_coverage`(조각 중 점이 덮은 비율), `puncta_overlap_area_um2`, `puncta_green_integrated`, `puncta_mean_area_um2`.
같은 세포 안 조각들은 서로 독립이 아니므로 조각 단위 p 값은 실제보다 작게(낙관적으로) 나옵니다.

---

## 3. 그림과 통계 읽는 법

- **Regression / correlation 히트맵**: 각 칸 = 두 지표 사이 Pearson r(또는 Spearman ρ). 빨강 = 양의 상관, 파랑 = 음의 상관. n < 3 이거나 값이 모두 같으면 계산하지 않음(–).
- **Groups: correlation**: 그룹 A, B 각각의 히트맵과 차이(Δ = B − A). `*` p < 0.05, `**` p < 0.01 (두 독립 상관계수의 Fisher z 검정).
  - 240개 칸을 동시에 검정하므로 p < 0.05 칸이 약 12개(5 %)는 우연으로도 나옵니다. 별 하나보다 **같은 방향으로 묶여 나오는 패턴**과 **이미지 단위에서도 유지되는지**를 보세요.
  - 세포 수가 적으면(n ≈ 15) r 값 하나는 ±0.5 정도 흔들릴 수 있습니다.
- **Groups: mean / median**: 지표별 그룹 비교. 점 = 세포(또는 이미지 평균), 막대 = 평균 ± SD / SEM / 95 % CI 또는 중앙값 ± IQR. 두 그룹은 Welch t와 Mann-Whitney U, 세 그룹 이상은 ANOVA와 Kruskal-Wallis.
  - *Unit = images* 로 바꾸면 이미지 하나를 한 번만 셉니다. 같은 이미지 세포들은 비슷하기 쉬우므로, 결론은 이미지 단위에서도 확인하는 것이 안전합니다.
- **Green+ / green−**: *Cells* 선택에서 green+ 만 고르면 발현 세포만 비교합니다. 발현 세포끼리 비교할 때 보통 이 설정을 씁니다.
