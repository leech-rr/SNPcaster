# SNPcaster
**English version is [here](#snpcasterenglish)**<br>
SNPcasterは次世代シーケンサーから得たショートリードデータに対し、以下の処理を行う解析パイプラインです。次の二つのプログラムから構成されています。
- Grape_qc_assembly
  - 品質チェック + アセンブリ
- SNPcaster
  - SNP解析 + 系統樹作成

## 講習会について
SNPcasterの使い方講習会を実施予定です。<br>
実施が決まりましたら、[Github Discussions](https://github.com/leech-rr/SNPcaster/discussions/8)にて告知します。

## インストール方法
以下の手順でSNPcasterを起動してください。<br>
詳細なインストール手順や使い方は[マニュアル](/doc/manual/2025-05-20_SNPcaster_インストール_操作マニュアル.pdf)に記載しています。ご不明点がある場合はそちらをご覧ください。

### ソースコードのダウンロード

#### gitを使った方法
gitを使用できる環境であれば、ターミナル上で以下のコマンドを使ってgit経由でソースコードをダウンロードできます。
```
git clone https://github.com/leech-rr/SNPcaster.git
```

#### zipを直接ダウンロードする方法
[ここをクリック](https://github.com/leech-rr/SNPcaster/archive/refs/heads/main.zip)するとダウンロードできます。</br>
ダウンロードしたら、お好きなフォルダに解凍してください。

### Dockerによる起動
SNPcasterはDockerによって環境構築を行います。<br>
[ソースコードのダウンロード](#ソースコードのダウンロード)を実施後、以下のコマンドをターミナルで実行してください。<br>
※cdに続くフォルダパスは実行環境に合わせて変更してください。
```
cd SNPcaster
docker compose up -d --build
```
初回起動は、1～数時間かかります。<br>
Dockerのインストール方法などの詳細は、[マニュアル](/doc/manual/2025-05-20_SNPcaster_インストール_操作マニュアル.pdf)に記載してあるのでご覧ください。

### Jupyter notebookへのアクセス
起動後、お好きなブラウザで[http://localhost:59829](http://localhost:59829)にアクセスすると、Jupyter notebookにアクセス可能です。

<div align="left">
  <img src="/doc/readme/images/jupyter_access.png" alt="Jupyterへのアクセス" style="width: 400px; border: 1px solid gray;">
</div>

[Create Project](http://localhost:59829/lab/tree/CreateProject_jp.ipynb)から新規のプロジェクトを作成して解析を行ってください。

## 引用
論文等で本プログラムを用いた方法を記載する際には、次の論文をご参照ください。
```
Lee K, Iguchi A, Uda K, Matsumura S, Miyairi I, Ishikura K, Ohnishi M, Seto J, Ishikawa K, Konishi N, Obata H, Furukawa I, Nagaoka H, Morinushi H, Hama N, Nomoto R, Nakajima H, Kariya H, Hamasaki M, Iyoda S. 2021. Whole-genome sequencing of Shiga toxin-producing Escherichia coli OX18 from a fatal hemolytic uremic syndrome case. Emerg Infect Dis 27:1509-1512.
```

## License
SNPcasterは[GPL v3.0](/COPYING)でライセンスされています。

## サポート
問題や質問がある場合は、[GitHub Issuesページ](https://github.com/leech-rr/SNPcaster/issues)を確認してください。

# SNPcaster(English)
SNPcaster is an analysis pipeline that performs the following processes on short-read data obtained from next-generation sequencers.The system is composed of the following two programs.
- Grape_qc_assembly
  - Quality check + Assembly
- SNPcaster
  - SNP analysis + Phylogenetic tree construction

## Installation
Please follow the steps below to start SNPcaster.<br>
Detailed installation instructions and usage are described in the [manual](/doc/manual/2025-05-20_SNPcaster_Installation_Operation_Manual.pdf). Please refer to it if you have any questions.

### Downloading the Source Code

#### Using Git
If you have a git environment, you can download the source code via git using the following command in the terminal:
```
git clone https://github.com/leech-rr/SNPcaster.git
```

#### Downloading the ZIP File Directly
You can download the ZIP file [here](https://github.com/leech-rr/SNPcaster/archive/refs/heads/main.zip).</br>
After downloading, please extract it to a folder of your choice.

### Starting with Docker
SNPcaster uses Docker for environment setup.<br>
After [downloading the source code](#downloading-the-source-code), execute the following commands in the terminal.<br>
※ Please adjust the folder path following `cd` according to your environment.
```
cd SNPcaster
docker compose up -d --build
```
The initial startup may take from 1 to several hours.<br>
For details on Docker installation and other information, please refer to the [manual](/doc/manual/2025-05-20_SNPcaster_Installation_Operation_Manual.pdf).

### Accessing Jupyter Notebook
Once started, you can access the Jupyter notebook by visiting [http://localhost:59829](http://localhost:59829) in your preferred browser.

<div align="left">
  <img src="/doc/readme/images/jupyter_access.png" alt="Accessing Jupyter" style="width: 400px; border: 1px solid gray;">
</div>

Create a new project from [Create Project](http://localhost:59829/lab/tree/CreateProject_jp.ipynb) to perform analysis.

## Citation
When describing methods using this program in papers or other publications, please refer to the following paper:
```
Lee K, Iguchi A, Uda K, Matsumura S, Miyairi I, Ishikura K, Ohnishi M, Seto J, Ishikawa K, Konishi N, Obata H, Furukawa I, Nagaoka H, Morinushi H, Hama N, Nomoto R, Nakajima H, Kariya H, Hamasaki M, Iyoda S. 2021. Whole-genome sequencing of Shiga toxin-producing Escherichia coli OX18 from a fatal hemolytic uremic syndrome case. Emerg Infect Dis 27:1509-1512.
```

## License
SNPcaster is licensed under [GPL v3.0](/COPYING).

## Support
For issues or questions, please check the [GitHub Issues page](https://github.com/leech-rr/SNPcaster/issues).
