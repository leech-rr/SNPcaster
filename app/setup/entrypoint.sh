#!/bin/bash
MOUNTED_DIR=${MY_MOUNTED_DIR:-/home/snpcaster/notebook/project}
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
chown -R "$EXEC_USER:$EXEC_GROUP" "$MOUNTED_DIR"

# コマンド実行
echo "Starting with UID: $HOST_UID, GID: $HOST_GID as USER: $EXEC_USER"
exec gosu "$EXEC_USER" "$@"
