## Project folder setup & link to GitHub

<aside>
<img src="/icons/book_gray.svg" alt="/icons/book_gray.svg" width="40px" />

> **Features:**
> 
> - How to organize your project folders in an intuitive way
> - How to sync your scripts to a folder where you cloned your GitHub repo
>     - Unfortunately, not all files can be backed up on GitHub due to storage limitations. The scripts should be sufficient to regenerate any data.
> 
> **Usage:**
> 
> - Part A: Add organization subfolders within your project folder
> - Part B: Rsync your scripts to your github repo folder semi-automatically
</aside>

### A. Project Folders

To keep myself organized, I structure each project folder as follows: 

Now, within the folder, there will (inevitably) be more sub-directories. This is because a project is divided in multiple analyses.

```bash
/scratch/mfafouti/Transfection_assessment_Tianze
├── Cell_type_annotation
│   ├── archive  # all file clutter goes here 
│   ├── data     # processed data that you generated (and wish to keep)
│   ├── logs     # all slurm logs can be directed here
│   ├── out
│   └── scripts
├── CR_count
│   ├── archive
│   ├── data
│   ├── logs
│   ├── out
│   └── scripts
├── CR_make_reference
│   ├── archive
│   ├── data
│   ├── logs
│   ├── out
│   └── scripts
├── **figure_data**   # all data required to just recreate a figure - can be a life saver
└── Scanpy_analysis
    ├── archive
    ├── data
    ├── logs
    ├── out
    └── scripts

31 directories, 0 files
```

The **`scripts` folder in each subdir is the most important one** and the only one I really want to keep backed up. I sometimes include some figures or text files. 

### Bash command: How to add all those subfolders to your analysis dirs at once

```bash
parent_dir="/path/to/parent"  # change this to your parent directory

for dir in "$parent_dir"/*/; do
    mkdir -p "$dir"/archive "$dir"/data "$dir"/logs "$dir"/out "$dir"/scripts
done
```

### B. GitHub & Rsync

The idea with GitHub is to keep all your scripts backed up after each day of work. Let’s make this as automated as possible, so it’s easy to stay consistent

I keep a Github folder on my `/scratch` where each sub-directory is a respository

```bash
/scratch/mfafouti/Github
├── Mommybrain-PPD # repo for project A 
│   ├── convert_to_h5seurat.R
│   ├── Figures
│   ├── LICENSE
│   ├── README.md
│   ├── seurat_env_narval.yml
│   ├── Slide_seq
│   └── Slide_tags
└── Transfection-assessment-scRNAseq # repo for project B
    ├── LICENSE
    ├── README.md
    └── test.py

6 directories, 7 files
```

However, I avoid editing this folder directly. Since my scripts are within the project folder, I set up rsync of the scripts folder in each project to the GitHub repo folder.

I wrote this script that I called `auto_rsync` and added it to `~/bin/` so it’s easy to execute from anywhere

### Rsync script

This takes as arguments your `project_dir` and `repo_dir` once you execute it: 

```bash
#!/bin/bash
# Hourly rsync of project scripts into GitHub repo
# Usage: ./auto_rsync.sh /path/to/project /path/to/repo

PROJECT_DIR="$1"
REPO_DIR="$2"

if [ -z "$PROJECT_DIR" ] || [ -z "$REPO_DIR" ]; then
    echo "Usage: $0 /path/to/project /path/to/repo"
    exit 1
fi

echo "Starting hourly rsync from $PROJECT_DIR to $REPO_DIR"

while true; do
    echo "🔄 $(date): Syncing scripts..."
    
    # Find all scripts folders and rsync
    find "$PROJECT_DIR" -type d -name scripts | while read -r scripts_path; do
        rel_path="${scripts_path#$PROJECT_DIR/}"
        dest="$REPO_DIR/$rel_path"
        mkdir -p "$dest"
        rsync -av --delete "$scripts_path/" "$dest/"
    done

    echo "✅ $(date): Sync complete. Sleeping 1 hour..."
    sleep 3600
done
```

Make it executable:

```bash
chmod +x ~/bin/auto_rsync
```

If you execute this in a new terminal window (from a login node) and allow it to run while you’re working, all your changes will be synced every hour.

Example of execution:

```bash
~/bin/auto_rsync /scratch/mfafouti/Transfection_assessment_Tianze /scratch/mfafouti/Github/Transfection-assessment-scRNAseq
```

## Set up automation to back up all of today’s scripts on your GitHub repo

You can add files and commit messages manually.. but there’s a more efficient way to do this. 

### Add a script called `git_backup` to `~/bin/` :

Note: I have gone through a process of setting up a GPG key to verify my commits, but this is not necessary. 

You will have to make a few changes to this script if you haven’t set this up

```bash
#!/bin/bash
# Interactive Git backup script
# Lists all repos under /scratch/mfafouti/Github/ and lets you choose one

GITHUB_DIR="/scratch/mfafouti/Github"

# Find all directories with a .git folder
repos=()
while IFS= read -r -d $'\0'; do
    repos+=("$(basename "$REPLY")")
done < <(find "$GITHUB_DIR" -maxdepth 1 -mindepth 1 -type d -exec test -d "{}/.git" \; -print0)

# Check if any repos found
if [ ${#repos[@]} -eq 0 ]; then
    echo "No Git repos found in $GITHUB_DIR"
    exit 1
fi

# Show menu
echo "Select a repo to back up:"
select repo in "${repos[@]}"; do
    if [ -n "$repo" ]; then
        echo "You selected: $repo"
        REPO_PATH="$GITHUB_DIR/$repo"
        break
    else
        echo "Invalid selection. Try again."
    fi
done

# Navigate to repo
cd "$REPO_PATH" || exit 1

# Show git status
git status

# Stage all changes
git add -A

# Show staged changes
echo "----------------------------------"
git status

# Setup GPG + SSH agent
export GPG_TTY=$(tty)
gpgconf --launch gpg-agent
ssh-add ~/.ssh/id_ed25519 2>/dev/null

# Ask for commit message
echo "Enter commit message (or leave empty for default timestamp):"
read commit_msg
if [ -z "$commit_msg" ]; then
    commit_msg="Auto-backup: $(date '+%Y-%m-%d %H:%M:%S')"
fi

# Commit with GPG signing
git commit -S -m "$commit_msg" || echo "⚠️ Nothing to commit."

# Push
git push origin main
```

Make it executable:

```bash
chmod +x git_backup
```

Whenever you want to submit to GitHub, you will just need to do:

```bash
git_backup
```

It will prompt you to 

1. Pick the repo you wish to back up
2. Write a commit message (optional, defaults to date)
3. Authorize the commit by entering your password (optional, but I like my GitHub to show that I was active)

Then it will automatically push to main.

## Daily workflow after you set this up

1. Make sure you have cloned your repo within a folder called `Github` in your scratch (no worries if the files are purged, they will be on your github either way)
    
    ```bash
    git clone https/link/to/your/repo
    ```
    
     Then, set up the SSH connection:
    
    ## Setting up SSH to your GitHub Repo (done individually for each repo)
    
    ```bash
    cd folder/to/your/repo
    ```
    
    Try: 
    
    ```bash
    git remote set-url origin [git@github.com](mailto:git@github.com):your-username/your-repo.git
    ```
    
    Now check again: 
    
    ```bash
    git remote -v
    ```
    
    The output should look like: 
    
    ```bash
    origin  git@github.com:neurogeek03/Mommybrain-PPD.git (fetch)
    origin  git@github.com:neurogeek03/Mommybrain-PPD.git (push)
    ```
    
    If you have already generated your SSH key, you can stop here
    
    1. Generate an SSH key (if you don't have one already)
        
        ```bash
        ssh-keygen -t ed25519 -C "your_email@example.com"
        ```
        
    2. Add the public key to GitHub
        1. Ask to see your public key
            
            ```bash
            cat ~/.ssh/id_ed25519.pub
            ```
            
        2. Copy the key
        3. Go to GitHub > Settings > SSH & GPG Keys
        4. Click: `Create a new Authentication Key`
        5. Paste what you copied earlier 
    3. Try: 
        
        ```bash
        ssh -T git@github.com
        ```
        
        You should get something like:
        
    
    You only need to do this once per repo.
    
2. Open a new terminal window and decide which project you will be rsyncing and to which repo:
    
    ```bash
    ~/bin/auto_rsync path/to/your/project/folder path/to/your/GitHub/repo/folder
    ```
    
    You should see:
    
    ```bash
    Starting hourly rsync from /scratch/mfafouti/Transfection_assessment_Tianze to /scratch/mfafouti/Github/Transfection-assessment-scRNAseq
    🔄 Tue 26 Aug 2025 04:01:03 PM EDT: Syncing scripts...
    sending incremental file list
    
    sent 65 bytes  received 12 bytes  154.00 bytes/sec
    total size is 0  speedup is 0.00
    sending incremental file list
    
    sent 63 bytes  received 12 bytes  150.00 bytes/sec
    total size is 0  speedup is 0.00
    sending incremental file list
    
    sent 87 bytes  received 12 bytes  198.00 bytes/sec
    total size is 0  speedup is 0.00
    sending incremental file list
    
    sent 64 bytes  received 12 bytes  152.00 bytes/sec
    total size is 0  speedup is 0.00
    sending incremental file list
    
    sent 64 bytes  received 12 bytes  152.00 bytes/sec
    total size is 0  speedup is 0.00
    ✅ Tue 26 Aug 2025 04:01:04 PM EDT: Sync complete. Sleeping 1 hour...
    ```
    
3. At the end of the day, run `git_backup` and follow the instructions!