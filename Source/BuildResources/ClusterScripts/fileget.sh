webdav_token=$(<webdav_token.txt)
webdav_url=$(<webdav_url.txt)
curl --user $webdav_token:nopass $webdav_url/$1 --remote-name
