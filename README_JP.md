# SNPcaster
**English version is [here](/README.md)**<br>
SNPcasterは次世代シーケンサーから得たショートリードデータに対し、以下の処理を行う解析パイプラインです。次の二つのプログラムから構成されています。
- Grape_qc_assembly
  - 品質チェック + アセンブリ
- SNPcaster
  - SNP解析 + 系統樹作成

## インストール方法
以下の手順でSNPcasterを起動してください。<br>
詳細なインストール手順や使い方は[マニュアル](/doc/manual/SNPcaster_manual_Japanese.pdf)に記載しています。ご不明点がある場合はそちらをご覧ください。

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
Dockerのインストール方法などの詳細は、[マニュアル](/doc/manual/SNPcaster_manual_Japanese.pdf)に記載してあるのでご覧ください。

### Jupyter notebookへのアクセス
起動後、お好きなブラウザで[http://localhost:59829](http://localhost:59829)にアクセスすると、Jupyter notebookにアクセス可能です。

<div align="left">
  <img src="/doc/readme/images/jupyter_access.png" alt="Jupyterへのアクセス" style="width: 400px; border: 1px solid gray;">
</div>

[Create Project](http://localhost:59829/lab/tree/CreateProject_jp.ipynb)から新規のプロジェクトを作成して解析を行ってください。

## 引用
論文等で本プログラムを用いた方法を記載する際には、GitHub page (https://github.com/leech-rr/SNPcaster)または次の論文をご参照ください。
```
Lee K, Iguchi A, Uda K, Matsumura S, Miyairi I, Ishikura K, Ohnishi M, Seto J, Ishikawa K, Konishi N, Obata H, Furukawa I, Nagaoka H, Morinushi H, Hama N, Nomoto R, Nakajima H, Kariya H, Hamasaki M, Iyoda S. 2021. Whole-genome sequencing of Shiga toxin-producing Escherichia coli OX18 from a fatal hemolytic uremic syndrome case. Emerg Infect Dis 27:1509-1512.
```
```
Lee K., Iguchi A., Terano C., Hataya H., Isobe J., Seto K., Ishijima N., Akeda Y., Ohnishi M., and Iyoda S. Combined usage of serodiagnosis and O antigen typing to isolate Shiga toxin-producing Escherichia coli O76:H7 from a hemolytic uremic syndrome case and genomic insights from the isolate. 2024. Microbiology spectrum 12:e0235523.
```

## リリース通知について
SNPcasterの新バージョンリリースの通知を受け取りたい方は、GitHubからリリース通知を受信できるように設定してください。<br>
設定方法は以下の通りです。
- GitHubアカウントを持っていなければ、作成してください(無料)。
  - [作成手順(GitHub公式)](https://docs.github.com/ja/get-started/start-your-journey/creating-an-account-on-github)
- [Githubにサインイン](https://github.com/login)します。
- [SNPcasterのページ](https://github.com/leech-rr/SNPcaster)にアクセスします。
- 画面右上にある「Watch」 > 「Custom」をクリックします。

<div align="left">
  <img src="/doc/readme/images/watch_github1.png" alt="Watchの設定" style="width: 400px; border: 1px solid gray;">
</div>

- 「Releases」を選択して、「Apply」をクリックします。

<div align="left">
  <img src="/doc/readme/images/watch_github2.png" alt="WatchでReleasesの選択" style="width: 400px; border: 1px solid gray;">
</div>

- これで、新バージョンがリリースされた際に、githubに登録したEメールアドレス宛てに通知が届くようになります。

### 通知が届かない場合
- 通知設定がうまくいかない場合は、[GitHub公式ドキュメント](https://docs.github.com/ja/account-and-profile/managing-subscriptions-and-notifications-on-github/setting-up-notifications/configuring-notifications#configuring-your-watch-settings-for-an-individual-repository)をご確認ください。

## License
SNPcasterは[GPL v3.0](/COPYING)でライセンスされています。

## サポート
問題や質問がある場合は、[GitHub Issuesページ](https://github.com/leech-rr/SNPcaster/issues)を確認してください。

## 免責事項
このソフトウェアはGPLv3ライセンスの下で提供されます。本ソフトウェアは無保証であり、使用によって生じたいかなる損害に対しても、作者は責任を負いません。
