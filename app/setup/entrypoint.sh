#!/bin/bash

USER_ID=${LOCAL_UID:-9001}
GROUP_ID=${LOCAL_GID:-9001}
USER_NAME=user

echo "Starting with UID : $USER_ID, GID: $GROUP_ID"
# UIDの確認：なかったら新しいユーザーを作成
if cat /etc/passwd | grep ${USER_ID}:${GROUP_ID} 1>/dev/null 2>&1 ; then
    echo "UID '${USER_ID}' already exists"
	USER_NAME=$(cat /etc/passwd | grep ${USER_ID}:${GROUP_ID} | awk -F ':' '{print $1}')
else
    useradd -u ${USER_ID} -o -m ${USER_NAME}
fi
# GIDの確認：なかったら新しいグループを作成、ある場合は既存グループに所属
if cat /etc/group | grep ${GROUP_ID} 1>/dev/null 2>&1 ; then
    echo "GID '${GROUP_ID}' already exists"
	GROUP_NAME=$(cat /etc/group | grep :${GROUP_ID}: | awk -F ':' '{print $1}')
	usermod -G ${GROUP_NAME} ${USER_NAME} 
else
    groupmod -g ${GROUP_ID} ${USER_NAME}
fi

SNPCASTER_HOME_DIR=/home/jovyan/snpcaster
rsync -av --update /usr/local/src/snpcaster/notebook "${SNPCASTER_HOME_DIR}/"

export HOME="${SNPCASTER_HOME_DIR}"
chown -R $LOCAL_UID:$LOCAL_GID "${SNPCASTER_HOME_DIR}"
chown -R $LOCAL_UID:$LOCAL_GID /usr/local

exec /usr/sbin/gosu ${USER_NAME} "$@"

