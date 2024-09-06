get run stats
-------------
use mymdc api

### intput
- `run_id`
- update frequency
- google sheet id

needs:
    - `credentials_google.json` (google)
    - `credentials_mdc.json` (mdc)

- authentication scripts
    - google (requires firefox) `authenticate_google.py`
        - need credentials file
    
    - myMDC (update token as required) `authenticate_mymdc.py`
        - need UID and secret (store in file) 
        - get from website

- build run stats sheet `update_run_stats.py`
    - call myMDC 
        - requires `experiment_id` which I can only get from a run 
    - possibly use EuXFEL tools
    - update every T seconds

- send to google spread sheet `update_google_sheet.py'

- would be nice to add progress on analysis 
    - but how to tell when a process is running or finished?


### Get UID and Secret from
- https://in.xfel.eu/metadata/oauth/applications/131

### Get token:
```
curl --request POST \
     --data grant_type=client_credentials \
     --data client_id=4qVGB9Bk1PQV5IsLcWYOZH4CVvpzVPkfaVP6LO11EWk \
     --data client_secret=... \
     https://in.xfel.eu/metadata/oauth/token
```
output: `{"access_token":"wOxr0O0nrxHBBkHsd2ZMZWedsH3WEvUj7LC_FYG6Vy0","token_type":"Bearer","expires_in":1480,"created_at":1725499637}`

get experiment id for 7076 = 1052
from the log book https://in.xfel.eu/metadata/proposals/1052#proposal-general

### api call for run stats
```
curl -X 'GET' \
  'https://in.xfel.eu/metadata/api/runs?experiment_id=1052' \
  -H 'accept: application/json; version=1' \
  -H 'X-USER-EMAIL: morganaj@unimelb.edu.au' \
  -H 'Authorization: Bearer _o7Lmdv-49SWwQNXr0qDJr6IidLs7-PD8dVohdIlKWs'
```


### get run info
- need experiment ids
    - need proposal id
        - need proposal number

save in json file so we don't need so many api calls
also only reload token when required
