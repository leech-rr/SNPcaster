# SNPcaster
**日本語版は[こちら](/README_JP.md)**<br>
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
When describing methods using this program in papers or other publications, please refer to GitHub page (https://github.com/leech-rr/SNPcaster) or the following papers:
```
Lee K, Iguchi A, Uda K, Matsumura S, Miyairi I, Ishikura K, Ohnishi M, Seto J, Ishikawa K, Konishi N, Obata H, Furukawa I, Nagaoka H, Morinushi H, Hama N, Nomoto R, Nakajima H, Kariya H, Hamasaki M, Iyoda S. 2021. Whole-genome sequencing of Shiga toxin-producing Escherichia coli OX18 from a fatal hemolytic uremic syndrome case. Emerg Infect Dis 27:1509-1512.
```
```
Lee K., Iguchi A., Terano C., Hataya H., Isobe J., Seto K., Ishijima N., Akeda Y., Ohnishi M., and Iyoda S. Combined usage of serodiagnosis and O antigen typing to isolate Shiga toxin-producing Escherichia coli O76:H7 from a hemolytic uremic syndrome case and genomic insights from the isolate. 2024. Microbiology spectrum 12:e0235523.
```


## Release Notifications
To receive notifications about new version releases of SNPcaster, configure your GitHub account to receive release notifications.  
Follow these steps:  
- If you don’t have a GitHub account, create one (it’s free).  
  - [Account Creation Guide (GitHub Official)](https://docs.github.com/en/get-started/start-your-journey/creating-an-account-on-github)  
- [Sign in to GitHub](https://github.com/login).  
- Visit the [SNPcaster page](https://github.com/leech-rr/SNPcaster).  
- Click **Watch** > **Custom** in the top-right corner of the page.  

<div align="left">  
  <img src="/doc/readme/images/watch_github1.png" alt="GitHub Watch button setup" style="width: 400px; border: 1px solid gray;">  
</div>

- Select **Releases** and click **Apply**.  

<div align="left">  
  <img src="/doc/readme/images/watch_github2.png" alt="Selecting Releases in GitHub Watch settings" style="width: 400px; border: 1px solid gray;">  
</div>  

- You will now receive notifications at the email address registered with your GitHub account whenever a new version is released.

### If You’re Not Getting Notifications
- If notifications aren’t arriving, check the [GitHub Official Documentation](https://docs.github.com/en/account-and-profile/managing-subscriptions-and-notifications-on-github/setting-up-notifications/configuring-notifications#configuring-your-watch-settings-for-an-individual-repository) for guidance.

## License
SNPcaster is licensed under [GPL v3.0](/COPYING).

## Support
For issues or questions, please check the [GitHub Issues page](https://github.com/leech-rr/SNPcaster/issues).
