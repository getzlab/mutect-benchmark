import requests
import json

# credential for datacommons
with open("/home/qingzhang/Downloads/credentials.json") as json_file: 
     credential = json.load(json_file)



## token for data commons
token = requests.post('https://nci-crdc.datacommons.io/user/credentials/api/access_token', 
    json=credential).json()
headers = {'Authorization': 'bearer '+ token['access_token']}


## token for GDC 
f = open("/home/qingzhang/Downloads/gdc-user-token.2019-10-24T14_36_35.270Z.txt", "r")
gdc_token = f.readlines()[0]
gdc_headers = {'X-Auth-Token': gdc_token}

## get UUID from manifest (need a better way to do this)
file_endpt = 'https://nci-crdc.datacommons.io/user/data/download/'


## from the GDC data commons, get access token
gdc_endpt = 'https://api.gdc.cancer.gov/legacy/data/'
file_uuid = '736a8e90-85ec-4007-b34a-1bf823eec6fc'
request_bai = "?&related_files=true"
gdc_response = requests.get(gdc_endpt+file_uuid+request_bai, 
    headers = gdc_headers)
print(json.dumps(gdc_response.json(), indent=2))


response = requests.get(file_endpt + file_uuid+request_bai, headers=headers)
print(json.dumps(response.json(), indent=2))



gatk GetSampleName -R refs/Homo_sapiens_assembly19.fasta -O sth.txt -I https://storage.googleapis.com/gdc-tcga-phs000178-controlled/KIRP/DNA/WXS/BCM/ILLUMINA/TCGA-UZ-A9PJ-01A-11D-A382-10_Illumina.bam?GoogleAccessId=zhangq07-3793@dcf-prod.iam.gserviceaccount.com&Expires=1571868331&Signature=kc1y5B5sg63opiiozoqVsnaY74xeT7apPnV71g%2BqOz2ZgvFS1uzrcmaf1hfcds6GMQgaU1ArTLNjYEA5s%2FEUl0Y1cbV3SSVcNZeSwqiKIiV8BwAHG34uPViP%2FG8xMxYiyFd5zA%2FAsK%2Fpm6SaIP2x3JVd%2F2%2FdAvXxa6rxgQD1nXJk8uIWFPn7adUMB%2BW3tUAFEFaoQU95RHifLkQzt0OwnherRECo%2FLSF0dMcLVoJ%2Fz5MrOuxM%2Fd%2BokcX6R%2BNzOPvCApM6gbTOXa5y24udiXxQo0IxatMqdlixD%2F3YEPib5Owhlr021BJH%2BmD%2FLE3FJ0cFTia8DKJ%2FBEmAW%2BQHNvZ3g==


# Data Download via API Endpoint Request:
durl = 'https://gen3.commons.io/api/v0/submission/<program>/<project>/export?format=bam&ids=' + ids[0:-1] # define the download url with the UUIDs of the records to download in "ids" list
dl = requests.get(durl, headers=headers)
print(dl.text) # display response

https://portal.gdc.cancer.gov/repository?
facetTab=files&
filters=%7B"op"%3A"and"%2C"content"%3A%5B%7B"op"%3A"in"%2C"content"%3A%7B
"field"%3A"cases.project.project_id"%2C
"value"%3A%5B"TCGA-KIRP"%5D%7D%7D%2C%7B"op"%3A"in"%2C"content"%3A%7B"field"%3A
"files.data_category"%2C"value"%3A%5B"Sequencing%20Reads"%5D%7D%7D%5D%7D&searchTableTab=files


&return_type=manifest
https://storage.googleapis.com/gdc-tcga-phs000178-controlled/KIRP/DNA/WXS/BCM/ILLUMINA/TCGA-UZ-A9PJ-01A-11D-A382-10_Illumina.bam?GoogleAccessId=zhangq07-3793@dcf-prod.iam.gserviceaccount.com&Expires=1571873978&Signature=l5R45egbczHXpxYB9lkN60miE3TTR2YIOQD5S4huIxn%2FiKWb4D%2FBpxFh36WjTYn0%2BgWNTUkdNZkzEyaU49CXAT2wHQj6nvY4xJuPtjU34M2ccDOzoK2nMp8K5ujJDDRDeEhcayjKicI2y8dXv7r5p4RqfEH2rRx5AZvjaTWoeHjBkjxef3u73ASc%2ByfHPtIRgzS6UhtYUfXkyIwhE34yGJSZlBGhmuwar9mtX914QgjGp0LSj2zFXFrGRwhqyEvm2j%2FUxlaqdhtZlA8cNmQt9jUiBZev2p%2B3MfJEYZ0PLMpy3bGhcQSaZssBh5k7UWjakj%2BaanIgH4RL4zsZwwh9Ew==

# https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/


files_endpt = 'https://api.gdc.cancer.gov/legacy/files'
params = {'fields':'cases.submitter_id,file_id,file_name'}
response = requests.get(files_endpt, params = params)


print(json.dumps(response.json(), indent=2))



query = {
    "op":"and",
    "content":[
        {
            "op":"in",
            "content":{
                "field":"files.data_category",
                "value":"TCGA-KIRP"
            }
        },
        {
            "op":"=",
            "content":{
                "field":"files.data_type",
                "value":"Aligned Reads"
            }
        }
    ]
}

files_endpt = 'https://api.gdc.cancer.gov/legacy/files'
params = {'filters':json.dumps(query),
    'fields':'cases.submitter_id,file_id,file_name'} 
    #'format':'tsv', 'return_type':'manifest'}
response = requests.get(files_endpt, params = params)
print(response)

print(json.dumps(response.json(), indent=2))



https://storage.googleapis.com/gdc-tcga-phs000178-controlled/KIRP/DNA/WXS/BCM/ILLUMINA/TCGA-UZ-A9PJ-01A-11D-A382-10_Illumina.bam?GoogleAccessId=zhangq07-3793@dcf-prod.iam.gserviceaccount.com&Expires=1571875572&Signature=l%2BCBY3dafI7MAfiYaRofhbe4qXDI6jXY1w6BAr4ZgPOFciwcMxCXiJyV98WsM1awjXXlfOrc2n%2BeEF0qGALfTHZikI4pSPFbX2MD1OxxjxH0obMwKh3g8P2Ux39qi%2BpZaS6U7d8A8al54LwfVAOqapSg8cFInNQLK8sWj92JNGVNFpAu5gTtW9yvAhmvsmdnpkfsrsT4QDM0AoZTmz9Lr%2FtN3IJUIGBqjCZXNEldfS%2F8Tk4NbA5mrrIsPb1aRN7760D5kWwASQ7%2FsqzrN07W3I4ajw1gSGEtdcWUAebwE2UggWb%2F7QoedrOXgh9%2Bx0kfVJzsWzC0acpdnJD4yHXutA==


# this 
https://storage.googleapis.com/gdc-tcga-phs000178-controlled/KIRP/DNA/WXS/BCM/ILLUMINA/TCGA-UZ-A9PJ-01A-11D-A382-10_Illumina.bam?GoogleAccessId=zhangq07-3793@dcf-prod.iam.gserviceaccount.com&Expires=1571876212&Signature=am1GBlZAGxwkh0qUHN7dlflGPY030mJLA2YxAxMydI9pJuPiJlMMRxVONPeKyFkRKnZAUDRXxxyNLt9%2F5Kt9No4Ax%2FJL6571EneAl0IOUwlr8bWqTferyOmFHWnJ44jDDMYtDS7VRIhe4dpVVqi9NOr3CYxuquqHQq4KBFTq9pYP5EsaJemptpGQCU5cGrQokwswrDb6l5YB9eI8H8JF9g6kQCZusLKR4Ltpof2eEym%2BIR0D9ckA%2FCGTLGBcF2%2BF73%2BKym6BvWg12CWCcLUUj%2FR%2BtdOn54I2Qk5RNQqPmYA2OH1box%2BUJkBEhplPXCAjC5poAjxfWra3Adiwy7iXmQ==




https://storage.googleapis.com/gdc-tcga-phs000178-controlled/KIRP/DNA/WXS/BCM/ILLUMINA/TCGA-UZ-A9PJ-01A-11D-A382-10_Illumina.bam?GoogleAccessId=zhangq07-3793@dcf-prod.iam.gserviceaccount.com&Expires=1571875183&Signature=LD6%2B%2Fp5OZzWwL3aXauoW0YFL9jdR6Q8KGCEJIJi3rSynmZ1eCl74S5lMzudfWhnipB3PvmqkkyIFnACv9bXGO8Da%2B8ABW0fZLykK2vzMx2xVzAxy0sJP5ZEz5QVaKKV%2FOzZVe8fIB8CJG3ypJ3bu3WmaMatMn%2BsDSXbqAQq8nzKVRxcXMub4wClFQZh%2B%2FCfnEmRIvK4xpR%2FGiJtZttdHt1YawuIZnMhRaiegCpht9bTvXtG%2BhBHt0e1tVFsauRoEFAmmkst6vK5OULdAtQZM5Wol%2BRrQdTRKuLA7hP5w7cluDnjm7Xm6OdQgVh6ljGNUE0b5Th0Dvp3ggoi%2FcnT5kg==