import subprocess
import uuid
# Generate a UUID4
 
import subprocess
import json
import uuid
import os
 
 
# Download the configuration JSON file from Google Cloud Storage
gcs_path = "gs://metro_tk_kplayground/control_json_can_128.json"
local_json_path = "/tmp/control_json.json"
 
print(f"Downloading configuration from {gcs_path}...")
subprocess.run(
    ["gcloud", "storage", "cp", gcs_path, local_json_path],
    check=True, capture_output=True, text=True
)
 
 
# Load the JSON configuration
with open(local_json_path, 'r') as f:
    config = json.load(f)
print(config)
dims = config.get("dims", "128x128x128")
add_radon = str(config.get("add_radon", False)).lower()
 
add_smooth = str(config.get("add_smooth", False)).lower()
additive_noise = str(config.get("additive_noise", 0.10))
 
is_can = str(config.get("is_can", False)).lower()
randomize = str(config.get("randomize", False)).lower()
variable_spacing = str(config.get("variable_spacing", False)).lower()
# Generate a unique ID
num_exec = int(config.get("num_exec", 1))
json_folders_path = config.get("json_folders_path", " ") # metro_tk_main/control_jsons
 
# Load filenames from Google Cloud Storage path
json_files = []
 
if(json_folders_path == ""):
    json_folders_path=" "
 
 
if (json_folders_path != " "):
    try:
        print(f"Listing files in {json_folders_path}...")
        result = subprocess.run(
            ["gcloud", "storage", "ls", json_folders_path],
            check=True, capture_output=True, text=True
        )
        # Parse the output to get file paths
        file_paths = result.stdout.strip().split('\n')
        # Filter out empty lines
        json_files = [path for path in file_paths if path]
        print(f"Found {len(json_files)} files in the GCS path")
    except subprocess.CalledProcessError as e:
        print(f"Failed to list files in GCS path: {e}")
        print(f"Error output: {e.stderr}")
   
    # Remove the GCS prefix from each file path
    json_files = [path.replace('gs://metro_tk_json/cans/', '') for path in json_files]
    # Filter out empty strings from json_files
    json_files = [file for file in json_files if file]
    num_exec=len(json_files)
    print(f"num_exec {num_exec}")    
print(f"json_folders_path {json_folders_path}   \n  json_files {json_files}\n")
 
 
for i in range(num_exec):
    try:
        # unique_id = str(uuid.uuid4())
        file = str(json_files[i])
        # unique_id = unique_id[12:-5]
        unique_id = dims + '_' + file[:-5]
 
       
        print("jjjjj is_can",is_can)
        # Choose the right Julia script based on type
        script_name = "in_docker_organized/main_create_phantom_ionic_chamber.jl"
        if config.get("is_can", False):
            script_name = f"in_docker_organized/main_create_phantom_can.jl"
           
       
        print(f"julia {script_name} {dims} {add_radon} {variable_spacing} {unique_id} {randomize} {add_smooth} {additive_noise}")
       
        # Execute the appropriate Julia script with parameters
        if(json_folders_path==" "):
            result = subprocess.run(
                ['julia', script_name, dims, add_radon, variable_spacing,unique_id,randomize,add_smooth,additive_noise],
                check=True, capture_output=True, text=True
            )
        else:
           
            # Download the JSON config file from Google Cloud Storage
            json_file_path = f"/tmp/{json_files[i]}"
            print(f"Downloading {json_files[i]} to {json_file_path}...")
            subprocess.run(
                ["gcloud", "storage", "cp", f"gs://metro_tk_json/cans/{json_files[i]}", json_file_path],
                check=True, capture_output=True, text=True
            )
            with open(json_file_path, 'r') as json_file:
                json_data = json.load(json_file)
            print(json_data)
 
            print("json_file_path:")
            print(json_file_path)
           
            result = subprocess.run(
                ['julia', script_name, dims, add_radon, variable_spacing,unique_id,randomize,add_smooth,additive_noise,json_file_path],
                check=True, capture_output=True, text=True
            )
        # Print the output
        # print(result.stdout)
       
    except subprocess.CalledProcessError as e:
        print(f"Command execution failed: {e}")
        print(f"Error output: {e.stderr}")
    except json.JSONDecodeError as e:
        print(f"Failed to parse JSON: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
    finally:
        # Clean up the temporary JSON file
        # if os.path.exists(local_json_path):
        #     os.remove(local_json_path)
        print("finally")
 
#.devcontainer/in_docker_organized/coordinate_phantom_create.py