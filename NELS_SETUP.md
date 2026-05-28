# NeLS Storage Setup

This project uses [NeLS storage](https://nels-docs.readthedocs.io/en/latest/index.html)
to share data files. The NeLS remote folder
`~/Projects/NMBU_threeD_Sandve_2024/Mucor/data` maps to the local `data/`
directory in this repository.

Data files are **not** downloaded wholesale. Instead, each analysis provides a
`download_data.sh` script that fetches only the files it needs (see
[Downloading data for an analysis](#4-downloading-data-for-an-analysis)).

---

## 1. Get your credentials from the NeLS portal

Log in to <https://nels.bioinfo.no/> and open **My Profile**.

**SSH key** — find the *SSH Private Key* section and click **Download Key**.
The file will be named `<username>.key` (e.g. `u15b2.key`). The portal also
shows the exact `chmod` command to run:

```
chmod 600 <username>.key
```

**SSH config** — find the *Example usage 2* section. It shows a ready-made
`~/.ssh/config` snippet (see step 3 below).

**Username** — your NeLS username is shown throughout that
page. Note it differs from your FEIDE/institutional login.

---

## 2. Install the SSH key

Move the downloaded key into `~/.ssh/` and set permissions:

```bash
mkdir -p ~/.ssh
mv ~/Downloads/<username>.key ~/.ssh/<username>.key
chmod 600 ~/.ssh/<username>.key
```

---

## 3. Add the SSH config

Open `~/.ssh/config` (create it if it does not exist) and paste the snippet
from the *Example usage 2* section on your **My Profile** page. It will look
like this (with your actual username substituted):

```sshconfig
## SSH gateway to nels-storage
Host nels.bastion
        HostName login.nels.elixir.no
        User <username>
        IdentityFile ~/.ssh/<username>.key
        ForwardAgent no
        IdentitiesOnly yes

## nels-storage host
Host nels.storage
        User <username>
        HostName data.nels.elixir.no
        IdentityFile ~/.ssh/<username>.key
        ProxyJump nels.bastion
        ForwardAgent no
        IdentitiesOnly yes
```

Set permissions on the config file:

```bash
chmod 600 ~/.ssh/config
```

Verify the connection:

```bash
ssh nels.storage
```

You should land in your NeLS home directory. Exit with `Ctrl-D`.

---

## 4. Downloading data for an analysis

Each analysis directory under `analysis/` contains a `download_data.sh` script
that lists the files it needs and delegates to the shared helper
`scripts/nels_download.sh`.

To download data for, say, the synteny analysis:

```bash
bash analysis/synteny/download_data.sh
```

Files are written to `data/<analysis>/...` relative to the repository root.
The `data/` directory is listed in `.gitignore` — data files are never
committed to the repository.

### Uploading results back to NeLS

After generating results you want to share, upload them with:

```bash
rsync -avh --progress data/<analysis>/ \
    nels.storage:Projects/NMBU_threeD_Sandve_2024/Mucor/data/<analysis>/
```

Add `-n` (`--dry-run`) to preview changes first.

