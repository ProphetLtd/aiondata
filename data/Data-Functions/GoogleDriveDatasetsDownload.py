from google.oauth2.credentials import Credentials
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request
import os
from googleapiclient.discovery import build  # Import build to create the service
from googleapiclient.http import MediaIoBaseDownload
import io


def authenticate_drive(TOKEN_PATH, CREDENTIALS_PATH, SCOPES):
    creds = None

    # Load token.json if it exists (and contains valid credentials)
    if os.path.exists(TOKEN_PATH):
        creds = Credentials.from_authorized_user_file(TOKEN_PATH, SCOPES)

    # If no valid credentials available, start the OAuth flow
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file(CREDENTIALS_PATH, SCOPES)
            creds = flow.run_local_server(port=0)
        # Save the new credentials for future use
        with open(TOKEN_PATH, "w") as token_file:
            token_file.write(creds.to_json())

    return creds


def list_files_in_folder(service, folder_id):
    query = f"'{folder_id}' in parents"
    results = (
        service.files()
        .list(
            q=query,
            pageSize=10,
            supportsTeamDrives=True,
            supportsAllDrives=True,
            includeItemsFromAllDrives=True,
            fields="nextPageToken, files(id, name)",
        )
        .execute()
    )
    items = results.get("files", [])
    if not items:
        print("No files found.")

    else:
        return items


def download_file(service, file_id, file_name, dir="."):
    request = service.files().get_media(fileId=file_id)
    fh = io.BytesIO()
    downloader = MediaIoBaseDownload(fh, request)
    done = False
    while not done:
        status, done = downloader.next_chunk()
        print(f"Download {int(status.progress() * 100)}%.")
    fh.seek(0)

    # Save to file
    with open(f"{dir}/{file_name}", "wb") as f:
        f.write(fh.read())
    print(f"File {file_name} downloaded successfully.")


if __name__ == "__main__":

    # Create a Project in the Google Cloud Console
    # Enable the Google Drive API
    # Create a "Dekstop App" OAuth 2.0 Client IDs in the Google Cloud Console
    # Allow Scopes: https://www.googleapis.com/auth/drive.readonly and other google drive scopes as needed
    # Download the credentials as a JSON file and save it in the auth folder
    # token.json will be created after the first authentication
    # The token.json file will store the user's credentials

    # Define the scopes required for accessing Google Drive
    SCOPES = ["https://www.googleapis.com/auth/drive.readonly"]
    # Paths to your OAuth credentials
    CREDENTIALS_PATH = (
        "auth/ProPhet_DownloadDrive_oAuth.json"  # Path to your web-style credentials
    )
    TOKEN_PATH = "auth/token.json"

    print("Authenticating...")

    # Authenticate and create the service
    creds = authenticate_drive(TOKEN_PATH, CREDENTIALS_PATH, SCOPES)
    service = build("drive", "v3", credentials=creds)  # Create the Drive API service

    print("Authenticated successfully.")

    # List and download the files in the Datasets folder
    folder_id = "1xFUzEv39GgPfjVlXYTD12rNyS81Kreks"  # Datasets folder
    files = list_files_in_folder(service, folder_id)
    dataset_dir = "local data/Datasets/"

    # download the files if they dont exist
    for file in files:
        if not os.path.exists(dataset_dir + file["name"]):
            download_file(service, file["id"], file["name"], dir=dataset_dir)
        else:
            print(f"File {file['name']} already exists.")
