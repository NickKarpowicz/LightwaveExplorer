webdav_token=$(<webdav_token.txt)
webdav_url=$(<webdav_url.txt)
curl --user $webdav_token:nopass $webdav_url/$1.zip --remote-name
unzip $1 -d $1
rm $1.zip
cd $1
sbatch $1.slurmScript

