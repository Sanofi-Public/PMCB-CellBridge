#!/bin/bash
while getopts a:t:i:d: flag
  do
    case "${flag}" in
      a) api=${OPTARG};;
      t) token=${OPTARG};;
      i) input=${OPTARG};;
      d) dest=${OPTARG};;
      *) echo "Invalid option: -$flag" ;;
  esac
done

if [ ! -f sb ]
then
    wget https://igor.sbgenomics.com/downloads/sb/linux-amd64/sb
fi
chmod +x sb

# echo "Seven Bridges API endpoint: $api"
# echo "Seven Bridges authorization token: $token"
# echo "input: $input"
# echo "destination: $dest"

cat <<EOF > credentials
[default]
api_endpoint   = $api
auth_token     = $token
advance_access = true
EOF

mkdir -p $HOME/.sevenbridges
mv credentials $HOME/.sevenbridges
./sb whoami
./sb upload start $input --destination "sanofi-research/$dest"
rm sb