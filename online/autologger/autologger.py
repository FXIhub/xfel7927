import update_run_stats
import update_google_sheet
import time

proposal_number = 7076
spreadsheet_id = "1cVhYDWk6-e_bODs7UhCC4EGg-iHNZApd31AjlywESao"
sheet_id       = '1278468409'

run_table_obj = update_run_stats.Run_table(proposal_number)

while True :
    # update run table
    headings, run_table = run_table_obj.update()
    
    # post to google spread sheet
    update_google_sheet.batch_update_values(spreadsheet_id, sheet_id, headings, run_table)
    
    print('sleeping')
    time.sleep(5)
