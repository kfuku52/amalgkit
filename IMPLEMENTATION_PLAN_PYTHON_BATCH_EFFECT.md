# Python化バッチ補正・正規化 実装計画

## 目的

AMALGKIT 全体から R 依存を段階的に排除し、最終的に runtime で `R` / `Rscript` / R package を不要にする。
本計画では以下を最終対象とする。

- `--batch_effect_alg sva` の Python 実装
- `--batch_effect_alg combatseq` の Python 実装
- `--batch_effect_alg ruvseq` の Python 実装
- `cstmm` で使っている `edgeR::calcNormFactors(method="TMM")` 相当の Python 実装
- `merge.r` の Python 置換
- `prepare_tables.r` / `wsfilter.r` / `finalize.r` の Python 置換
- `cross_species_filter.r` の Python 置換
- exclusion plot を含む PDF/可視化出力の Python 置換
- R 向け runtime probe、CI、テスト、設定受け渡しの削除
- 追加機能として `--batch_effect_alg latent_glm` を実装し、非負値保証の補正経路を提供する

本計画は「実装しながら更新する運用文書」として扱う。フェーズ完了時、仕様変更時、受け入れ条件変更時に必ず更新する。

## Final Outcome (2026-03-08)

- repo 全体の runtime / CLI / CI から R 依存を削除済み
- `merge`, `cstmm`, `wsfilter`, `csfilter`, `finalize`, `cross_species_filter` は Python 実装で動作する
- batch correction の `sva`, `combatseq`, `ruvseq` は Python backend で動作する
- plot / PDF 出力は `matplotlib` ベースへ移行済み
- `check_rscript()`、R package probe、R smoke job、`.r` エントリーポイント、R reference helper は削除済み
- repo 内の tracked `.r` / `.R` ファイルは削除済み
- `python -m pytest -q` は `912 passed`、`python -m compileall amalgkit tests` も通過済み
- 今後の作業は R-free 化ではなく、Python 実装の精度改善や新機能追加として扱う

## 最重要要件

- 最終到達点は「AMALGKIT 全体の R フリー」である
- `sva` は既存 R backend と高い一致性を確認する
- `edgeR` 必要部分は既存 R 実装と高い一致性を確認する
- `ruvseq` / `combatseq` も最終的には Python 化する
- `latent_glm` は parity 対象ではなく、新機能として非負値保証と完走性を重視する
- 既存 CLI と既存出力ディレクトリ構造をなるべく壊さない
- Python 側の描画は `matplotlib` を標準とし、plot/PDF の最終実装は `matplotlib` ベースに統一する

## 非目標

- `RUVSeq` / `ComBat_seq` / `sva` の full package port
- byte-for-byte 一致の保証
- 初期段階での全 `.r` ファイル一括削除

## 現在の主要差し込み点

- `amalgkit/finalize.r`
  - `batch_effect_subtraction()`
  - `run_batch_effect_step()`
- `amalgkit/per_species_tables.py`
  - `run_per_species_r_script()`
  - `config_map` の引き渡し
- `amalgkit/cli_parser.py`
  - `finalize` サブコマンドの batch-effect 関連オプション
- `amalgkit/cstmm.py`
  - `cstmm.r` 呼び出し
- `amalgkit/cstmm.r`
  - `edgeR::calcNormFactors(method="TMM")`
- `amalgkit/merge.r`
- `amalgkit/prepare_tables.r`
- `amalgkit/wsfilter.r`
- `amalgkit/cross_species_filter.r`
- `amalgkit/filter_utils.py`
  - `save_exclusion_plot_pdf()`
- `amalgkit/runtime_utils.py`
  - `check_rscript()`
- `amalgkit/cli_utils.py`
  - R / Rscript / R package probe
- `.github/workflows/tests.yml`
  - R smoke/integration job

## 現在の前提整理

- `combatseq` は現在 `sva::ComBat_seq()` を使っているため、runtime から `sva` package を完全に消すには `combatseq` も Python 化する必要がある
- `ruvseq` は現在 `edgeR` と `RUVSeq` の両方に依存している
- `merge`, `prepare_tables`, `wsfilter`, `finalize`, `cross_species_filter` の可視化と集計ロジックもまだ R に残っている
- したがって「`sva` を Python 化」「`edgeR` 必要部分を移植」「`ruvseq` / `combatseq` を Python 化」に加えて、「各 `.r` エントリーポイント自体を Python に置き換える」必要がある
- ただし作業順は一括ではなく、`sva` → `combatseq` → `edgeR/TMM` → `ruvseq` → per-species Python 化 → cross-species Python 化 → plot Python 化 → runtime/CI cleanup → `latent_glm` が安全

## 運用ルール

- 各フェーズ開始時に本ファイルの `Status` を更新する
- 仕様が変わったら `Decision Log` を更新する
- parity 条件が緩和・強化されたら `受け入れ条件` を更新する
- 新しいテスト fixture を追加したら `テスト資産` に記録する
- ブロッカーが出たら `Open Questions` に追記する

## Status

- Phase 0: completed
- Phase 1: completed
- Phase 2: completed
- Phase 3: completed
- Phase 4: completed
- Phase 5: completed
- Phase 6: completed
- Phase 7: completed
- Phase 8: completed
- Phase 9: completed
- Phase 10: completed
- Phase 11: completed
- Phase 12: completed

### Progress Snapshot

- R-free migration target is complete
- repo 全体で Python-only runtime が成立し、`.r` / `.R` source は削除済み
- full pytest suite は `919 passed`
- CI の R smoke job は削除済み
- `merge`, `cstmm`, `prepare_tables`, `wsfilter`, `finalize`, `cross_species_filter` の Python path は production に接続済み
- exclusion plot, tau histogram, merge PDF 群, batch compare PDF は `matplotlib` 実装へ移行済み
- `check_rscript()`、`temporary_r_config`、`r_config.py`、R package probe は削除済み
- `latent_glm` backend は production path に接続済みで、manual/auto の backend test、per-species integration、CLI smoke を通過済み
- 以降の残タスクは新機能追加や数値精度改善であり、R 依存除去と `latent_glm` 実装は完了している

- CLI / config 伝搬は完了
- `finalize.r` から Python SVA backend を呼ぶ bridge は実装済み
- `sva_backend=python` の explicit `n.sv=0` と auto-estimated `n.sv=0` は end-to-end で通る
- `num.sv(method='be')` 相当の初期 Python 実装を追加済み
- `irwsva.build()` 相当の初期 Python 実装を追加済み
- `num.sv(method='leek')` fallback の初期実装を追加済み
- `sva_backend=python` の manual positive `n.sv` は `finalize.r` 経由で完走する
- `sva_backend=python` 経路では `sva` R package を不要化済み
- direct backend parity test を追加し、raw synthetic matrix では R reference と corrected matrix 一致を確認済み
- transformed FPKM/log2p1 duplicate-case で `num.sv(be)` unresolved 時に `leek` fallback する制御を修正済み
- `finalize.r` end-to-end parity は manual positive `n.sv` と auto-estimation の両方で R backend 一致を確認済み
- SVA 関連の targeted suite は `36 passed` で通過済み
- Python `combatseq` backend を `inmoose.pycombat_seq` wrapper として追加済み
- `combatseq_backend=python` の direct backend test と `finalize.r` end-to-end parity を balanced fixture で確認済み
- `edgeR::calcNormFactors(method=\"TMM\")` 相当の Python コアを `normalization_tmm.py` として追加済み
- TMM round1 auto-reference 選定、median reference 再計算、multiple median refs を含む direct parity を `edgeR` と確認済み
- `cstmm` の Python backend は single-species に加えて multi-species まで接続済み
- `cstmm_backend=python` は synthetic fixture と R backend parity を通過済み
- `ruvseq_backend=python` の初期版を `batch_effect_ruvseq.py` と `finalize.r` bridge に接続済み
- `RUVr` 数式部分は `RUVSeq::RUVr` との direct parity を確認済み
- `finalize.r` の `ruvseq_backend=python` end-to-end test を追加済み
- `ruvseq` 前段の p-value ranking / residual 生成は OLS fallback から `statsmodels` GLM 近似へ引き上げ済み
- `ruvseq` backend 全体の R oracle helper を追加し、manual/auto の `resolved_ruv_k` 一致を test で固定済み
- auto-k の score tie-break は `k=0` との同点時に最小の正の `k` を優先するよう安定化済み
- `seqUQ` 行列は edgeR 近似から `EDASeq::betweenLaneNormalization(which=\"upper\")` 相当に切り替え済み
- residual / p-value 生成は Poisson-first + effective library size に更新済みで、balanced fixture では R と一致
- `finalize.r` の balanced `ruvseq` fixture で R/Python parity を確認済み
- single-group `sample_group` の `ruvseq_design_failed` no-op も Python/R parity を test で固定済み
- `save_exclusion_plot_pdf()` は `matplotlib` 実装へ置換済みで、Rscript なしで PDF を出力できる
- `merge.r` 相当の plot 生成は `merge_plots.py` へ移植済みで、`merge.py` から Python 実装を直接呼ぶ
- `merge_main()` は `check_rscript()` なしで動作し、`merge` command から Rscript 必須を外した
- `tests/test_merge.py` に merge PDF 出力の Python regression を追加済み
- `per_species_common.py` に `sample_group_mean` / `sample_group_to_tau` の parity helper を追加し、R reference test を通過済み
- `per_species_common.py` に round summary / final curation summary writer を追加し、unit test で固定済み
- `per_species_outputs.py` に index付きTSV writer、correlation statistics、tau histogram PDF の Python helper を追加済み
- correlation statistics は R reference と parity test を通過済みで、tau histogram PDF は `matplotlib` 実装で regression test を追加済み
- `batch_effect_common.py` に `initialize_batch_info` / metadata annotation / batch-effect summary writer を追加し、unit test で固定済み
- `per_species_finalize_python.py` を追加し、`finalize.r` 相当の `skip_curation` / `disable_auto_outlier_filter` 経路を Python 実装へ移植済み
- `per_species_tables.py` は supported な `finalize` species job について Python worker を優先し、Rscript check を不要化済み
- `finalize.py` は current default path と同じく `disable_auto_outlier_filter=True` を明示し、Python worker selection と整合済み
- `tests/test_per_species_finalize_python.py` と `tests/test_filter_commands.py` に、Rscript を禁止した状態での Python-only finalize end-to-end test を追加済み
- 現在の targeted regression suite は `187 passed` で通過済み
- 未完了項目は SVA fixture 拡充、`--sva_backend python` / `--combatseq_backend python` / `--ruvseq_backend python` の default 化判断、`ruvseq` の control selection / GLM residual の R 近似精度向上、`merge.r` 以外の残存 plot/PDF の Python 化、per-species / cross-species の本体 Python 化、runtime/CI 全体からの R probe 削除

## 全体アーキテクチャ方針

### 基本方針

- 初期段階では `finalize.r` をオーケストレータとして残し、batch correction の中核計算から Python 化する
- 中間段階では per-species / cross-species 処理本体を Python モジュールへ移し、R は比較用 backend に縮退させる
- 最終段階では PDF、summary、metadata 出力も含めて Python 化し、`.r` エントリーポイントを削除する
- Python backend は TSV/JSON ベースの明示的な入出力を持つ
- backend ごとの差分は共通 schema に吸収する
- Python の描画実装は `matplotlib` を中核にし、`seaborn` などの高水準ラッパーには依存しない

### 新規モジュール候補

- `amalgkit/batch_effect_common.py`
- `amalgkit/batch_effect_io.py`
- `amalgkit/batch_effect_sva.py`
- `amalgkit/batch_effect_combatseq.py`
- `amalgkit/batch_effect_ruvseq.py`
- `amalgkit/batch_effect_latent_glm.py`
- `amalgkit/normalization_tmm.py`
- `amalgkit/per_species_pipeline.py`
- `amalgkit/cross_species_pipeline.py`
- `amalgkit/plotting.py`
- `amalgkit/exclusion_plot.py`
- `amalgkit/batch_effect_cli.py`
- `amalgkit/scripts/run_batch_backend.py`

## 新CLI仕様

### 既存維持

- `--batch_effect_alg no|sva|ruvseq|combatseq`

### 追加

- `--sva_backend r|python`
- `--batch_effect_alg latent_glm`
- `--latent_family poisson|nb`
- `--latent_k auto|INT`
- `--latent_k_max INT`
- `--latent_max_iter INT`
- `--latent_tol FLOAT`

### デフォルト方針

- 初期段階では `--sva_backend r` を default とする
- parity が安定したら `--sva_backend python` を default に切り替える
- `latent_glm` は opt-in の新機能として導入する

## 共通I/O設計

### Python backend 入力

- counts matrix TSV
- metadata TSV
- options JSON

### Python backend 出力

- corrected matrix TSV
- backend summary JSON
- optional latent/surrogate matrix TSV
- optional diagnostics TSV

### summary JSON 共通キー

- `backend`
- `method`
- `skip_reason`
- `stable`
- `corrected_runs`
- `uncorrected_runs`
- `resolved_sva_nsv`
- `resolved_sva_B`
- `resolved_ruv_k`
- `resolved_ruv_controls`
- `negative_values_before_clip`
- `negative_values_after_clip`

backend ごとに不要キーは `null` を許容する。

## フェーズ別計画

### Phase 0: Oracle 固定とテスト基盤

#### 目的

- 現行 R 実装を参照系として凍結する
- 以後の Python 実装比較の基準を作る

#### 作業

- `tests/fixtures/batch_effect/` を新設
- `sva`, `combatseq`, `ruvseq`, `cstmm` 用の小さな fixture を作る
- `merge`, `wsfilter`, `finalize`, `csfilter` 用の出力 fixture も最小限作る
- 以下のケースを最低限含める
  - confounded design
  - single sample
  - singleton batch
  - `n.sv=0`
  - auto `B`
  - manual `B`
  - zero-heavy matrix
  - negative clipping 対象
  - sample_group が 1 群しかないケース
- 既存 PDF 出力の存在、ファイル名、ページ数、メタデータ列の oracle を定める
- 現行 R backend の出力を snapshot として保存する

#### 追加テスト候補

- `tests/test_finalize_r_oracle.py`
- `tests/test_cstmm_r_oracle.py`

#### 受け入れ条件

- すべての fixture について、現行 R 実装の出力を再生成できる
- 参照値の生成手順が README か本計画に記録される

### Phase 1: CLI と設定伝搬

#### 目的

- 新 backend を安全に選べるようにする

#### 変更対象

- `amalgkit/cli_parser.py`
- `amalgkit/per_species_tables.py`
- `amalgkit/finalize.r`
- `tests/test_cli_help.py`
- `tests/test_curate.py`

#### 作業

- `finalize` の CLI に新オプション追加
- `config_map` に `sva_backend` と `latent_*` を追加
- `finalize.r` の config 読み込み部に対応項目を追加
- help text を更新

#### 受け入れ条件

- 既存コマンドのデフォルト挙動が変わらない
- 新オプションが `-h` に表示される

### Phase 2: Python SVA core 実装

#### 目的

- `--batch_effect_alg sva --sva_backend python` を成立させる

#### 変更対象

- 新規 `amalgkit/batch_effect_common.py`
- 新規 `amalgkit/batch_effect_io.py`
- 新規 `amalgkit/batch_effect_sva.py`
- 新規 `amalgkit/scripts/run_batch_backend.py`
- 必要に応じて `amalgkit/subprocess_utils.py`

#### 実装対象

- design matrix 構築補助
- `resolve_sva_B_value`
- `estimate_sva_nsv_at_B`
- `resolve_sva_parameters`
- `sva()` 相当の surrogate variable 推定
- `cleanY()` 相当の residualization
- skip reason と summary 情報の生成

#### 重要仕様

- 参照実装は Bioconductor `sva()` / `num.sv()` / `cleanY`
- `method='be'` を優先し、失敗時に `leek` へ fallback
- `resolved_sva_nsv` と `resolved_sva_B` の決定ロジックは現行 `finalize.r` と一致させる
- `num.sv(method='be')` が unresolved を返すケースでは `leek` に fallback し、transformed duplicated-column case でも R と同じ `auto_leek` 解決を再現する

#### 数値実装方針

- runtime 依存は原則 `numpy` + `scipy`
- `statsmodels` は検証用に限る
- 線形代数のランク、投影、残差計算は明示的に実装し、BLAS 差の影響を小さくする

#### 受け入れ条件

- fixture 上で `resolved_sva_nsv` が R と一致
- fixture 上で `resolved_sva_B` が R と一致
- corrected matrix が tight tolerance で一致
- `skip_reason` が R と一致

### Phase 3: `finalize.r` への Python SVA bridge 統合

#### 目的

- 既存 finalize 出力を壊さず Python SVA を利用可能にする

#### 変更対象

- `amalgkit/finalize.r`
- 新規 `amalgkit/batch_effect_r_bridge.R` は必要なら追加
- `tests/test_finalize_integration.py`

#### 作業

- `batch_effect_subtraction()` に `sva_backend == "python"` 分岐追加
- Python helper を `system2()` で呼び出す
- Python 出力を既存の `sva1` 互換 list に再構成
- 既存の `draw_sva_summary()` と `.RData` 保存との整合を取る

#### 受け入れ条件

- 既存の PDF 生成と summary 出力が壊れない
- `--sva_backend r` と `--sva_backend python` を並存運用できる

### Phase 4: SVA parity 検証の強化

#### 目的

- SVA の R 置換を実務上受け入れ可能なレベルまで固める

#### 追加テスト

- `tests/test_batch_effect_sva.py`
- `tests/test_finalize_parity.py`

#### 比較対象

- `resolved_sva_nsv`
- `resolved_sva_B`
- `skip_reason`
- corrected matrix
- metadata の `batch_corrected`
- metadata の `batch_alg_used`
- summary TSV の主要列

#### 比較方法

- corrected matrix: `numpy.testing.assert_allclose`
- SV matrix: 符号反転や回転不定性があるため subspace 比較
- TSV/metadata: 列ごとの exact compare を基本とする

#### 受け入れ条件

- CI 上で parity test が安定して通る
- 差分が出るケースは既知差分として文書化される

### Phase 5: Python Combat-seq 実装

#### 目的

- `--batch_effect_alg combatseq` を Python backend に移す

#### 実装候補

- 既存 Python port を採用して vendor/freeze する
- もしくは repo 内実装を追加する

#### 変更対象

- 新規 `amalgkit/batch_effect_combatseq.py`
- `amalgkit/scripts/run_batch_backend.py`
- `amalgkit/finalize.r`
- `tests/test_batch_effect_combatseq.py`

#### 仕様

- singleton batch は現行と同様に uncorrected のまま保持
- `group` 指定ありで失敗したら group なし fallback
- corrected/unorrected run の扱いは現行 summary と整合させる

#### 受け入れ条件

- singleton batch あり/なしで既存挙動を維持
- corrected matrix が R と十分近い
- `skip_reason` と corrected_runs が期待通り

### Phase 6: edgeR 必要部分の Python 実装

#### 目的

- `cstmm` と `ruvseq` 前段のために edgeR 必要部分を置き換える

#### 実装範囲

- TMM normalization factor
- lib size の取り扱い
- upper-quartile normalization
- 必要最小限の CPM 相当
- RUVSeq 前段で必要な GLM 準備のうち edgeR 相当部分

#### 変更対象

- 新規 `amalgkit/normalization_tmm.py`
- 必要に応じて `amalgkit/edge_like_glm.py`
- `amalgkit/cstmm.py`
- 新規 `amalgkit/cstmm_python.py`
- `tests/test_tmm_edger_parity.py`
- `tests/test_cstmm.py`

#### parity 対象

- `calcNormFactors(method="TMM", refColumn=NULL)`
- median ref column 選定
- `calcNormFactors(..., refColumn=median_index)`
- normalization factors
- metadata に書かれる `tmm_library_size` / `tmm_normalization_factor`

#### 受け入れ条件

- fixture 上で TMM factor が R と tight tolerance で一致
- `cstmm` 出力 counts と metadata が R と一致

#### 現状

- TMM factor parity は `edgeR` と一致確認済み
- `cstmm` の single-species / multi-species Python backend は実装済み
- multi-species は synthetic fixture と R backend parity を通過済み
- `ruvseq` 前段で使う upper-quartile normalization / CPM / residual 準備の初期 Python 実装は `batch_effect_ruvseq.py` に着手済み

### Phase 7: Python RUVSeq 実装

#### 目的

- `--batch_effect_alg ruvseq` を Python backend に移す

#### 実装範囲

- 現 repo が使っている `RUVr` 経路に限定
- control gene auto selection
- upper-quartile normalization
- GLM fitting
- deviance residuals
- `RUVr` 相当の unwanted factor 推定
- `k` auto selection

#### 変更対象

- 新規 `amalgkit/batch_effect_ruvseq.py`
- `amalgkit/batch_effect_runner.py`
- `amalgkit/finalize.r`
- `amalgkit/finalize.py`
- `amalgkit/per_species_tables.py`
- `amalgkit/cli_parser.py`
- `tests/test_batch_effect_ruvseq.py`

#### 注意点

- これは `RUVSeq` package 全移植ではなく repo 使用範囲の移植
- `controls` の選び方、`k` 選択、スコア算出は現行 `finalize.r` を参照する

#### 受け入れ条件

- `resolved_ruv_k` が期待レンジで再現する
- corrected matrix が R と十分近い
- `ruvseq_k_zero`、`ruvseq_fit_failed` などの skip reason が再現される

#### 現状

- `ruvseq_backend=python` の CLI/config/`finalize.r` bridge は実装済み
- `RUVr` 本体は R source に沿って Python 化し、direct parity test を追加済み
- control selection / auto-k / residual 生成は Python 近似実装を追加済み
- residual 生成と control ranking は `statsmodels` GLM 近似を優先し、失敗時のみ軽量 fallback を使う構成に更新済み
- backend 全体を比較する R oracle helper と pytest を追加済み
- auto-k は R と同じ `resolved_ruv_k` を返すよう tie-break を安定化済み
- `seqUQ` の補正行列は `EDASeq` upper normalization 相当へ置換済みで、manual/auto の corrected matrix も R に近づいている
- balanced fixture では backend 直接比較・`finalize.r` 統合比較の両方で corrected matrix parity を通過済み
- `finalize.r` の Python `ruvseq` backend は manual `k=1` fixture で end-to-end 完走済み
- まだ一般 fixture 全体での parity は未達で、今後の主作業は dispersion 推定をさらに寄せて適用範囲を広げること

### Phase 8: per-species R pipeline の Python 置換

#### 目的

- `prepare_tables.r` / `wsfilter.r` / `finalize.r` の実処理を Python 化する

#### 変更対象

- 新規 `amalgkit/per_species_pipeline.py`
- `amalgkit/per_species_tables.py`
- `amalgkit/wsfilter.py`
- `amalgkit/finalize.py`
- 必要に応じて `amalgkit/metadata_utils.py`
- `tests/test_curate.py`
- `tests/test_filter_commands.py`
- `tests/test_finalize_integration.py`

#### 実装範囲

- transform logic
- sample-group mean
- tau
- mapping rate filter
- within-species outlier filter
- metadata update
- curation summary
- batch-effect summary

#### 受け入れ条件

- `prepare_tables.r`, `wsfilter.r`, `finalize.r` を Python 経路で置換できる
- 主要 TSV 出力が既存 R 版と互換である
- `per_species_tables.py` から R スクリプト呼び出しが不要になる

### Phase 9: cross-species R pipeline の Python 置換

#### 目的

- `cross_species_filter.r` の本体を Python 化する

#### 変更対象

- 新規 `amalgkit/cross_species_pipeline.py`
- `amalgkit/cross_species_filter.py`
- `amalgkit/csfilter.py`
- `tests/test_csca.py`
- `tests/test_filter_commands.py`

#### 実装範囲

- ortholog table 集約
- averaged / unaveraged expression assembly
- missing value handling
- PCA / MDS / t-SNE 用データ生成
- cross-species outlier flagging
- metadata 出力

#### 受け入れ条件

- `cross_species_filter.py` から R スクリプト呼び出しが不要になる
- `csfilter` の metadata と excluded.tsv が既存版と互換である

### Phase 10: plot / PDF の Python 置換

#### 目的

- すべての可視化出力の R 依存をなくす

#### 変更対象

- 新規 `amalgkit/plotting.py`
- 新規 `amalgkit/exclusion_plot.py`
- `amalgkit/merge.py`
- `amalgkit/filter_utils.py`
- `amalgkit/finalize.py`
- `amalgkit/wsfilter.py`
- `amalgkit/csfilter.py`
- `tests/test_merge.py`
- `tests/test_filter_commands.py`

#### 実装範囲

- exclusion plot
- merge plot
- cstmm plot
- per-species plot
- cross-species plot

#### 描画方針

- すべての主要 plot/PDF は `matplotlib` で実装する
- heatmap、scatter、boxplot、histogram、dendrogram、legend layout も `matplotlib` ベースで統一する
- 必要に応じて `scipy` の補助関数は使うが、描画ライブラリの中核は `matplotlib` とする
- `seaborn` は必須依存にしない

#### 受け入れ条件

- `save_exclusion_plot_pdf()` が Python 実装になる
- `merge.r` の出力 PDF が Python から生成される
- R なしで全主要 PDF が生成される

### Phase 11: runtime / test / CI から R 依存を除去

#### 目的

- runtime probe、テスト、CI、配布設定から R 依存を取り除く

#### 変更対象

- `amalgkit/runtime_utils.py`
- `amalgkit/cli_utils.py`
- `.github/workflows/tests.yml`
- `tests/test_curate_integration.py`
- `tests/test_finalize_integration.py`
- `tests/r_smoke/`
- `setup.py`
- `README.md`

#### 作業

- `check_rscript()` を不要化または削除
- runtime banner から R / Rscript / R packages を除外
- R smoke job を削除
- R 前提 integration test を Python backend 前提へ置換
- ドキュメントから R 前提を削除

#### 受け入れ条件

- CI が R なしで通る
- runtime で `Rscript` を探しに行かない
- README と help が新状態に一致する

### Phase 12: R エントリーポイント削除と最終クリーンアップ

#### 目的

- リポジトリ全体を R フリーにする

#### 変更対象

- `amalgkit/*.r`
- `amalgkit/r_config.py`
- R bridge があればそれも削除
- 関連テストとドキュメント

#### 作業

- 不要になった `.r` ファイル削除
- `temporary_r_config()` と DCF 受け渡し削除
- package manifest の整理
- dead code 削除

#### 最終受け入れ条件

- runtime で `R`, `Rscript`, R package が不要
- リポジトリに機能実行に必要な `.r` ファイルが残っていない
- 主要 CLI が Python 実装のみで完走する

### Phase 13: latent_glm 実装

#### 目的

- 非負値を保証する新しい batch correction backend を提供する

#### 状態

- completed

#### 変更対象

- 新規 `amalgkit/batch_effect_latent_glm.py`
- `amalgkit/per_species_finalize_python.py`
- `amalgkit/batch_effect_runner.py`
- `amalgkit/batch_effect_common.py`
- `amalgkit/cli_parser.py`
- `tests/test_batch_effect_latent_glm.py`
- `tests/test_finalize_integration.py`
- `tests/test_per_species_finalize_python.py`
- `tests/test_batch_effect_runner.py`

#### 実装方針

- raw count を入力
- Poisson/NB log-link
- `counts ~ group + latent + offset`
- latent factor と係数を交互最適化
- 出力は adjusted counts

#### 受け入れ条件

- 最終出力まで完走する
- 補正後 matrix に負値を出さない
- synthetic data で batch signal が減り、group signal が一定以上保持される

#### 完了メモ

- backend 本体は `amalgkit/batch_effect_latent_glm.py` に実装済み
- `--batch_effect_alg latent_glm` は per-species finalize worker と runner に接続済み
- `--latent_family`, `--latent_k`, `--latent_k_max`, `--latent_max_iter`, `--latent_tol` は CLI から使用可能
- auto-k は residual SVD energy を使う heuristic で決定する
- corrected counts は常に非負へ clip され、summary TSV に `resolved_latent_k`, `latent_family`, `latent_iterations`, `latent_objective`, `latent_converged` を出力する
- `python -m pytest -q` は `919 passed`
- `Rscript` を PATH から外した `finalize --batch_effect_alg latent_glm` CLI smoke は通過済み

## テスト資産

今後追加する想定のテスト。

- `tests/test_batch_effect_sva.py`
- `tests/test_finalize_parity.py`
- `tests/test_batch_effect_combatseq.py`
- `tests/test_tmm_edger_parity.py`
- `tests/test_batch_effect_ruvseq.py`
- `tests/test_batch_effect_latent_glm.py`
- `tests/test_per_species_pipeline.py`
- `tests/test_cross_species_pipeline.py`
- `tests/test_plotting_outputs.py`

fixture 候補。

- `tests/fixtures/batch_effect/sva_confounded/`
- `tests/fixtures/batch_effect/sva_nsv_zero/`
- `tests/fixtures/batch_effect/combatseq_singleton/`
- `tests/fixtures/batch_effect/ruvseq_auto_k/`
- `tests/fixtures/cstmm/multispecies_tmm/`

## CI計画

- 初期段階では R backend と Python backend の両方を CI で回す
- parity test は小さい fixture のみ CI 常時実行
- 重い比較は nightly または手動ジョブ化を検討
- `sva_backend=python` が安定したら default 切り替え
- plot/per-species/cross-species 置換後に R job を削除する
- 最終段階で `cli_utils.py` の runtime banner から R / Rscript / R packages を削除する

## リスクと対策

### 数値差

- リスク: BLAS/LAPACK 差、乱数差、SVD 符号反転
- 対策: matrix 自体と subspace を分けて検証する

### SVA auto 推定差

- リスク: permutation 実装差で `n.sv` が揺れる
- 対策: fixture を固定し、seed と permutation 仕様を明文化する

### RUVSeq の再現性

- リスク: `edgeR` 相当の GLM 再現が難しい
- 対策: まず repo 使用範囲に限定し、control selection と residuals を重点的に比較する

### Combat-seq 依存

- リスク: 外部 Python port の挙動差や保守性
- 対策: vendor 固定または wrapper を薄く保ち swap 可能にする

### 既存 PDF 出力との整合

- リスク: `sva_out` 構造が微妙に違って plot が壊れる
- 対策: `finalize.r` 側で既存互換 list を組み直す

### 全体 R フリー化の終盤での回帰

- リスク: `.r` 削除時に見落とした entrypoint が残る
- 対策: `rg` による `Rscript`, `.r`, `check_rscript`, `temporary_r_config` のゼロ件確認を受け入れ条件にする

## Decision Log

- 初期段階では `finalize.r` 全廃はしない
- `sva` は Python backend を主軸にし、R backend は parity 用に残す
- `combatseq` も Python 化対象に含める
- Python `combatseq` backend の初期実装には `inmoose.pycombat_seq` を採用する
- Python TMM 実装は `edgeR::calcNormFactors.default(method=\"TMM\")` の数式と reference 選定を直接移植する
- `ruvseq` は Python 化対象に含めるが、repo 使用範囲の `RUVr` 相当に限定する
- 最終ゴールは batch-effect 中核の Python 化ではなく repo 全体の R フリー化とする
- Python 側の plot/PDF 実装は `matplotlib` に統一する
- `latent_glm` は parity 対象ではなく新機能として扱う
- transformed duplicated-column case の `be -> leek` fallback は Python 側で再現済み

## Open Questions

- `combatseq` は既存 Python port を採用するか、自前実装するか
- `ruvseq` の GLM 部分をどこまで edgeR 互換に寄せるか
- `sva` parity の許容誤差を fixture ごとに変えるか
- `.RData` の保存形式をどこまで既存互換にするか
- どの時点で `--sva_backend python` を default に切り替えるか
- cross-species plot の完全互換をどこまで要求するか

## 作業開始順

最初のマイルストーンは以下。

1. Phase 0 の fixture 固定
2. Phase 1 の CLI/config 追加
3. Phase 2 の Python SVA core 実装
4. Phase 3 の `finalize.r` bridge 接続
5. Phase 4 の parity test 合格

この後のマイルストーンは以下。

6. Phase 5-7 で `combatseq`、`edgeR/TMM`、`ruvseq` を Python 化
7. Phase 8-10 で per-species / cross-species / plot を Python 化
8. Phase 11-12 で runtime / CI / repository から R を除去
9. Phase 13 で `latent_glm` を追加し、全体を整理する
