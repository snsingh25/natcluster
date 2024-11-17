
# How To : Guide on how to use git

### 1. Clone the Main Repository
```
git clone <main-repo-url>
cd <main-repo-directory>
```

### 2. Add a Submodule
```
git submodule add <submodule-repo-url> <path/to/submodule-directory>
git commit -m "Add <submodule-directory> as a submodule"
```

### 3. Initialize and Update Submodules (After Cloning)
```
git submodule update --init --recursive
```

### 4. Pull Updates for Submodules
```
git submodule update --remote --merge
git commit -am "Update submodules to latest commit"
```

### 5. Clone a Repository with Submodules
```
git clone --recurse-submodules <main-repo-url>
```

### 6. If Submodules Aren't Cloned with the Main Repository
```
git submodule update --init --recursive
```
