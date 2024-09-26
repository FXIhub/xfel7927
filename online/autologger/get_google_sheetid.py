import authenticate_google
from googleapiclient.discovery import build
from googleapiclient.errors import HttpError

# p7076
#spreadsheetId = "1cVhYDWk6-e_bODs7UhCC4EGg-iHNZApd31AjlywESao"


# p7927
spreadsheetId = "15nfDVdlXoTenMnGrbUlyC005MD8xmpgZX0R5y1fnQeI"
creds = authenticate_google.get()

service = build("sheets", "v4", credentials=creds)

fields = 'sheets.properties'

result = (
    service.spreadsheets()
    .get(spreadsheetId = spreadsheetId, fields=fields)
    .execute()
)

print(result)
