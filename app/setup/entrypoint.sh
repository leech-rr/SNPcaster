#!/bin/bash
HOME_DIR=${MY_MOUNTED_DIR:-/home/snpcaster}
MOUNTED_DIR=${MY_MOUNTED_DIR:-${HOME_DIR}/notebook/project}
echo "MOUNTED_DIR: $MOUNTED_DIR"

# ホストのnotebooksのUIDとGIDを取得
HOST_UID=$(stat -c %u "${MOUNTED_DIR}")
HOST_GID=$(stat -c %g "${MOUNTED_DIR}")

# 既存ユーザーを検索
EXISTING_USER=$(getent passwd "$HOST_UID" | cut -d: -f1)
echo "EXISTING_USER: $EXISTING_USER"
if [ -n "$EXISTING_USER" ]; then
    # 既存ユーザーがいればそのユーザー
    EXEC_USER="$EXISTING_USER"
else
    # いなければsnpcasterのUIDを変更
    usermod -u "$HOST_UID" snpcaster || true
    EXEC_USER="snpcaster"
fi

# 既存グループを検索
EXISTING_GROUP=$(getent group "$HOST_GID" | cut -d: -f1)
echo "EXISTING_GROUP: $EXISTING_GROUP"
if [ -n "$EXISTING_GROUP" ]; then
    # 既存グループがあればそのグループ
    EXEC_GROUP="$EXISTING_GROUP"
else
    # いなければsnpcasterのGIDを変更
    groupmod -g "$HOST_GID" snpcaster || true
    EXEC_GROUP="snpcaster"
fi

# ユーザーのプライマリグループを設定
usermod -g "$EXEC_GROUP" "$EXEC_USER"

# パーミッションを調整
chown -R "$EXEC_USER:$EXEC_GROUP" "$HOME_DIR"

# --- custom_setup.sh の初回実行 ---
CUSTOM_SETUP_SCRIPT="/usr/local/bin/custom_setup.sh"
FLAG_FILE="/var/run/.custom_setup_done"
if [ -f "$CUSTOM_SETUP_SCRIPT" ]; then
    if [ ! -f "$FLAG_FILE" ]; then
        echo "[Entrypoint] Running custom_setup.sh for the first time..."
        chmod +x "$CUSTOM_SETUP_SCRIPT"
        bash "$CUSTOM_SETUP_SCRIPT"
        touch "$FLAG_FILE"
        echo "[Entrypoint] custom_setup.sh finished."
    else
        echo "[Entrypoint] custom_setup.sh has already been executed. Skipping."
    fi
fi
# -----------------------------------

# コマンド実行
echo "Starting with UID: $HOST_UID, GID: $HOST_GID as USER: $EXEC_USER"
exec gosu "$EXEC_USER" "$@"