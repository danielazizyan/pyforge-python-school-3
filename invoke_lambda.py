import boto3
import json

client = boto3.client('lambda')

event = {
    "names": ["Alice", "Bob"]
}

response = client.invoke(
    FunctionName='HelloStudentFunction',
    InvocationType='RequestResponse',
    Payload=json.dumps(event)
)

response_payload = json.loads(response['Payload'].read())

print(response_payload)
