import update_run_stats
import update_google_sheet

proposal_number = 7076
sheet_id = "1cVhYDWk6-e_bODs7UhCC4EGg-iHNZApd31AjlywESao"

while True :
    # update run table
    headings, run_table = update_run_stats.update_run_table(proposal_number)

    # post to google spread sheet
    update_google_sheet.batch_update_values(sheet_id, headings, run_table)
