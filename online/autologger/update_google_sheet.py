import authenticate_google
from googleapiclient.discovery import build
from googleapiclient.errors import HttpError
import json
import numpy as np


# help: https://developers.google.com/sheets/api/samples/writing
# https://developers.google.com/sheets/api/guides/values#python
# write data
def batch_update_values(spreadsheet_id, sheet_id, headings, values, col_start = 0, row_start = 0, value_input_option = "USER_ENTERED" ):
    """
    Creates the batch_update the user has access to.
    Load pre-authorized user credentials from the environment.
    TODO(developer) - See https://developers.google.com/identity
    for guides on implementing OAuth2 for the application.

    writes a list of lists values to sheets starting at column col_start and row row_start
    """
    creds = authenticate_google.get()
    
    data = [headings] + values
    
    # make input for batch update
    # I wish I could use service.spreadsheets().values().batchUpdate
    # but I can't figure out how to add a sheet id to this
    # https://stackoverflow.com/questions/64099894/python-google-sheets-api-make-values-update-and-sheet-properties-update-with-a-s
    rows = []
    for r in data:
        col = []
        for c in r:
            col.append({"userEnteredValue": ({"numberValue": c} if str(c).replace('.', '', 1).isdigit() else {"stringValue": c})})
            #col.append({"userEnteredValue": ({"stringValue": str(c)})})
        rows.append({"values": col})
    
    # https://stackoverflow.com/questions/50916422/python-typeerror-object-of-type-int64-is-not-json-serializable
    #rows = json.dumps(rows, cls=NpEncoder)

    body = {
    "requests": [
        {
            "updateCells": {
                "start": {
                    "sheetId": sheet_id,
                    "rowIndex": row_start,
                    "columnIndex": col_start
                },
                "rows": rows,
                "fields": "userEnteredValue"
            }
        }
    ]
    }
    
    
    #try:
    if True:
        service = build("sheets", "v4", credentials=creds)
        
        result = (
            service.spreadsheets()
            .batchUpdate(spreadsheetId=spreadsheet_id, body=body)
            .execute()
        )
        
        print(f"{(result.get('totalUpdatedCells'))} cells updated.")
        return result
    
    #except HttpError as error:
    #    print(f"An error occurred: {error}")
    #    return error

