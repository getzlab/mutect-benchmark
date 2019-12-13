

import requests
import json
with open("/home/qingzhang/Downloads/credentials.json") as json_file: 
         credential = json.load(json_file)
token = requests.post('https://nci-crdc.datacommons.io/user/credentials/api/access_token', 
        json=credential).json()
headers = {'Authorization': 'bearer '+ token['access_token']}

file_endpt = 'https://nci-crdc.datacommons.io/user/data/download/'

file_uuid = "3ae0c34c-d80d-430e-8cdc-0d30e85a8d47"
response = requests.get(file_endpt + file_uuid, headers=headers)
# get signed url here:
print(json.dumps(response.json(), indent=2))

# for the bams, we need to index by ourselfs - 
# or we can download directly from GDC repository 
# by curl
