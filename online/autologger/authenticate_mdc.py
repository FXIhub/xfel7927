from metadata_client.metadata_client import MetadataClient
import json
import os

# load credentials
credentials = json.load(open('credentials_mdc.json'))

def get(fnam = 'credentials_mdc.json'):
    credentials = json.load(open(fnam))
    
    client_id     = credentials['client_id']
    client_secret = credentials['client_secret']
    user_email    = credentials['email']
    token_url     = 'https://in.xfel.eu/metadata/oauth/token'
    refresh_url   = 'https://in.xfel.eu/metadata/oauth/token'
    auth_url      = 'https://in.xfel.eu/metadata/oauth/authorize'
    scope         = '' 
    base_api_url  = 'https://in.xfel.eu/metadata/api'
    
    proposal_number = 7076
    
    client_comm = MetadataClient(client_id = client_id,
                                 client_secret = client_secret,
                                 user_email = user_email,
                                 token_url = token_url,
                                 refresh_url = refresh_url,
                                 auth_url = auth_url,
                                 scope = scope,
                                 base_api_url = base_api_url)
                                 
    #all_proposal_runs = MetadataClient.get_proposal_runs(client_comm, proposal_number)
    #print(json.dumps(all_proposal_runs, indent=2))
    return client_comm
