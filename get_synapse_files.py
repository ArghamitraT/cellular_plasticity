#,'syn23626795'
file_list = [ 'syn23627528',
             'syn23627309',
             'syn23626843',
             'syn23626844',
             'syn25421678',
             'syn23627482',
             'syn23626823',
             'syn23626825',
            'syn23627552']

import synapseclient
syn = synapseclient.Synapse()
syn.login(authToken='eyJ0eXAiOiJKV1QiLCJraWQiOiJXN05OOldMSlQ6SjVSSzpMN1RMOlQ3TDc6M1ZYNjpKRU9VOjY0NFI6VTNJWDo1S1oyOjdaQ0s6RlBUSCIsImFsZyI6IlJTMjU2In0.eyJhY2Nlc3MiOnsic2NvcGUiOlsidmlldyIsImRvd25sb2FkIiwibW9kaWZ5Il0sIm9pZGNfY2xhaW1zIjp7fX0sInRva2VuX3R5cGUiOiJQRVJTT05BTF9BQ0NFU1NfVE9LRU4iLCJpc3MiOiJodHRwczovL3JlcG8tcHJvZC5wcm9kLnNhZ2ViYXNlLm9yZy9hdXRoL3YxIiwiYXVkIjoiMCIsIm5iZiI6MTcwMTExMTIyMywiaWF0IjoxNzAxMTExMjIzLCJqdGkiOiI0NDMwIiwic3ViIjoiMzQ4NzUwOCJ9.hY7AcFvYcP7Gcp4LeRJUtr1TV2tdVXA_tJoNSd6YaUz-WmroxlBRs5oLpT8sBp5mpWhVkllcm6lKUyKD6Hg4XqIsG24IMEIpRRaa9ePsu4EF9av0ShCuR049LAz_waRsfzDEW1Aa4xQcQoqsblPmg9_YeE00NWFDSto8a8ydX34iUIDF2LHoMJL-DFcKQEx-ZUXXeg0NZBk70xNttqLIvUsDAae1IT7vrGDe4cP3kJwQHdx6Do8CJfd6fJhkcagoAQDaRHGnOC9Rd6OPvGos9O8s5NunrA6hJNhnC-BjA1h2uLO1KaAaRv84BCLsfZ2gB_f460LThEjSGkyw_su83A', rememberMe=True)

for file in file_list:
    entity = syn.get(file,
    downloadLocation="/Users/arghamitratalukder/Library/CloudStorage/GoogleDrive-at3836@columbia.edu/My Drive/technical_work/Lineage_plasticity/data/MSK_tumor_data")
    print("done ", file)



