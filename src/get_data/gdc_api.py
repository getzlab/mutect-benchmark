

import requests
import json
with open("/home/qingzhang/Downloads/credentials.json") as json_file: 
         credential = json.load(json_file)
         #query = {'query':"""{project(first:0){project_id TCGA-COAD}}"""}
         #credential.update(query)
token = requests.post('https://nci-crdc.datacommons.io/user/credentials/api/access_token', 
        json=credential).json()
headers = {'Authorization': 'bearer '+ token['access_token']}

file_endpt = 'https://nci-crdc.datacommons.io/user/data/download/'
         # file_endpt = 'https://api.gdc.cancer.gov/legacy/data/'
#file_uuid = '736a8e90-85ec-4007-b34a-1bf823eec6fc'
file_uuid = "3ae0c34c-d80d-430e-8cdc-0d30e85a8d47"
response = requests.get(file_endpt + file_uuid, headers=headers)
print(json.dumps(response.json(), indent=2))

