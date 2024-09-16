import authenticate_google
from googleapiclient.discovery import build
from googleapiclient.errors import HttpError

spreadsheetId = "1cVhYDWk6-e_bODs7UhCC4EGg-iHNZApd31AjlywESao"
creds = authenticate_google.get()

service = build("sheets", "v4", credentials=creds)

fields = 'sheets.properties'

result = (
    service.spreadsheets()
    .get(spreadsheetId = spreadsheetId, fields=fields)
    .execute()
)

print(result)
