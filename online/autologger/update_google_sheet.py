import authenticate_google
from googleapiclient.discovery import build
from googleapiclient.errors import HttpError

# help: https://developers.google.com/sheets/api/samples/writing
# https://developers.google.com/sheets/api/guides/values#python
# write data
def batch_update_values(spreadsheet_id, headings, values, col_start = 0, row_start = 0, value_input_option = "USER_ENTERED" ):
    """
    Creates the batch_update the user has access to.
    Load pre-authorized user credentials from the environment.
    TODO(developer) - See https://developers.google.com/identity
    for guides on implementing OAuth2 for the application.

    writes a list of lists values to sheets starting at column col_start and row row_start
    """
    creds = authenticate_google.get()

    # make google range string
    assert(col_start < 26)
    col_start_string = chr(ord('@')+1+col_start)
    col_stop_string  = chr(ord('@')+1+col_start+len(values[0]))
    
    row_start_string = 2+row_start
    row_stop_string  = 2+len(values)
    
    range_name = f'{col_start_string}{row_start_string}:{col_stop_string}{row_stop_string}'
    

    col_start_string = chr(ord('@')+1+col_start)
    col_stop_string  = chr(ord('@')+1+col_start+len(headings))
    
    row_start_string = 1+row_start
    row_stop_string  = 1+row_start
    
    range_name_headings = f'{col_start_string}{row_start_string}:{col_stop_string}{row_stop_string}'
    
    try:
        service = build("sheets", "v4", credentials=creds)
        
        # perhaps values needs to be a list of lists
        data = [{'range': range_name_headings, 'values': [headings]},
                {'range': range_name, "values": values}]
        #data = [{'range': range_name, "values": values}]
         
        body = {"valueInputOption": value_input_option, "data": data}
        
        result = (
            service.spreadsheets()
            .values()
            .batchUpdate(spreadsheetId=spreadsheet_id, body=body)
            .execute()
        )
        
        print(f"{(result.get('totalUpdatedCells'))} cells updated.")
        return result
    
    except HttpError as error:
        print(f"An error occurred: {error}")
        return error
