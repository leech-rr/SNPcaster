#!/bin/bash -e

source .env

echo "" >> ~/.bashrc
echo "# Porxy setteing by set_proxy.sh" >> ~/.bashrc
echo "export http_proxy=\"${HTTP_PROXY}\"" >> ~/.bashrc
echo "export HTTP_PROXY=\"${HTTP_PROXY}\"" >> ~/.bashrc
echo "export https_proxy=\"${HTTPS_PROXY}\"" >> ~/.bashrc
echo "export HTTPS_PROXY=\"${HTTPS_PROXY}\"" >> ~/.bashrc
echo "export no_proxy=\"localhost, 127.0.0.0/8, ::1\"" >> ~/.bashrc
echo "export NO_PROXY=\"\$no_proxy\"" >> ~/.bashrc


echo "# Porxy setteing by set_proxy.sh" | sudo tee /etc/sudoers.d/proxy
echo "Defaults env_keep += \"http_proxy https_proxy no_proxy HTTP_PROXY HTTPS_PROXY NO_PROXY\"" | sudo tee -a /etc/sudoers.d/proxy
echo "Done!"
