import json
import os

# load credentials
credentials = json.load(open('credentials_mdc.json'))

# how to check if existing token is still valid?
# just load every time

def get():
    # get token
    token = os.popen(
    f"""
    curl --request POST \
         --data grant_type=client_credentials \
         --data client_id={credentials['client_id']} \
         --data client_secret={credentials['client_secret']} \
         https://in.xfel.eu/metadata/oauth/token
    """).read()
    
    creds = {'token': json.loads(token)['access_token'], 'email': credentials['email']}
    return creds
