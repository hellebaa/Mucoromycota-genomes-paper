# NeLS Storage Setup

This project uses [NeLS storage](https://nels-docs.readthedocs.io/en/latest/index.html)
to share data files. Raw and intermediate data live in the NeLS folder
`~/Projects/NMBU_threeD_Sandve_2024/Mucor` and are synchronized with the local
`data/` directory in this repository.

## 1. Get your credentials from the NeLS portal

Log in to <https://nels.bioinfo.no/> and open **My Profile**. From there,
download:

- your **SSH private key** (e.g. `nels-userkey-XXXX.txt`)
- the suggested **SSH config snippet** (shown on the same page)
- your NeLS **username** (note that it is *not* your FEIDE ID)

NeLS uses a jump host (`login.nels.elixir.no`) in front of the storage host
(`data.nels.elixir.no`), so two hops are required for every connection.

## 2. Install the SSH key

Move the key into your `~/.ssh/` directory and lock down the permissions:

```bash
mkdir -p ~/.ssh
mv ~/Downloads/nels-userkey-*.txt ~/.ssh/nels_key
chmod 600 ~/.ssh/nels_key
```

## 3. Add the SSH config

Open `~/.ssh/config` (create it if it does not exist) and paste the snippet
from the NeLS portal. It should look roughly like this — replace
`<username>` with your NeLS username:

```sshconfig
Host nels-login
    HostName login.nels.elixir.no
    User <username>
    IdentityFile ~/.ssh/nels_key
    IdentitiesOnly yes

Host nels
    HostName data.nels.elixir.no
    User <username>
    IdentityFile ~/.ssh/nels_key
    IdentitiesOnly yes
    ProxyJump nels-login
```

Set the permissions on the config file:

```bash
chmod 600 ~/.ssh/config
```

Verify the connection:

```bash
ssh nels
```

You should land in your NeLS home directory (`/nels/users/<username>`). Exit
with `Ctrl-D`.

## 4. Synchronize `data/` with NeLS

The shared data lives at `~/Projects/NMBU_threeD_Sandve_2024/Mucor` on NeLS
and mirrors to `data/` in this repository.

Make sure the local target exists:

```bash
mkdir -p data
```

### Download (NeLS → local)

```bash
rsync -avh --progress nels:Projects/NMBU_threeD_Sandve_2024/Mucor/ data/
```

### Upload (local → NeLS)

```bash
rsync -avh --progress data/ nels:Projects/NMBU_threeD_Sandve_2024/Mucor/
```

The trailing slashes matter: with them, the *contents* of the source folder
are copied into the destination folder.

Add `-n` (`--dry-run`) to preview changes before running for real, and
`--delete` if you want the destination to exactly mirror the source (removes
files that no longer exist on the source side — use with care).

### Interactive browsing

For one-off file operations you can also use `sftp`:

```bash
sftp nels
sftp> cd Projects/NMBU_threeD_Sandve_2024/Mucor
sftp> ls
```

Graphical clients such as FileZilla or Cyberduck work too — point them at
`data.nels.elixir.no`, set the jump host to `login.nels.elixir.no`, and use
the same private key.
