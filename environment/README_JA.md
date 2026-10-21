# RとStanの準備

モデル推定にはR、**CmdStan**、**C++コンパイラ**が必要です。

## 1. Rパッケージを用意する

作業ディレクトリをリポジトリ最上位にして実行します。

```r
source("environment/setup.R")
install_cantab_packages()
```

不足するパッケージをCRANとStan公式のR-universeからインストールします。既に入っているパッケージを自動で別の版に置き換えません。A′のモデルには`xbeta`を提供する版の`brms`が必要です。

`brms`が`xbeta`に対応していない場合は、新しいRセッションで更新します。

```r
install.packages("brms", repos = "https://cloud.r-project.org")
```

更新後はRを再起動し、準備と環境確認を再実行します。上記はCRANの現行版を導入する操作で、下記の検証環境を復元するものではありません。

## 2. CmdStanを用意する

未導入なら、OSに合うコンパイラを準備し、次を実行します。

```r
source("environment/setup.R")
install_cantab_packages(install_stan = TRUE)
```

コンパイラはWindowsならRの版に対応するRtools、macOSならC++のコマンドライン開発ツールを使います。実行時のツールチェーン確認と[CmdStanR公式手順](https://mc-stan.org/cmdstanr/articles/cmdstanr.html)の案内に従ってください。

既存のCmdStanを独自の場所に入れている場合は、`bin`フォルダを含むインストール先を指定します。

```r
stan_dir <- normalizePath("path/to/cmdstan", winslash = "/", mustWork = TRUE)
Sys.setenv(CMDSTAN = stan_dir)
cmdstanr::set_cmdstan_path(stan_dir)
```

環境変数`CMDSTAN`は、解析が起動する子Rプロセスにインストール先を引き継ぐために必要です。

## 3. 環境を確認する

```r
source("environment/check_environment.R")
```

検証に使用した環境はR 4.6.1、brms 2.23.0、mgcv 1.9-4、cmdstanr 0.9.0、posterior 1.7.0、CmdStan 2.39.0です。[パッケージ一覧](recorded_R_packages.csv)は参照用で、環境を自動復元するファイルではありません。版が違う場合は表を確認し、実行中・再開時は同じ環境を使います。

実行確認はLinuxで行っています。Windows/macOS実機の実行確認は含みません。環境を用意したら、[実行ガイド](../docs/USAGE_JA.md)の`mode = "check"`から始めてください。計算の条件は[設定ガイド](../docs/SETTINGS_JA.md)、ログや確認結果の場所は[診断ガイド](../docs/DIAGNOSTICS_JA.md)にあります。
