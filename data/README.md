# Data directory

Research data are intentionally excluded from version control. This prevents accidental redistribution of controlled, licensed, very large, or independently maintained resources.

Recommended local layout:

```text
data/
├── raw/
│   ├── peaks/
│   ├── expression/
│   ├── epigenomics/
│   ├── genetics/
│   └── reference/
├── interim/
└── processed/
```

Set the environment variable `TCF4_DATA_DIR` to use an external data location, or copy `config/example_config.yml` to `config/config.yml` and edit the paths. Do not commit participant-level, controlled-access, or licensed source data.

The current scripts originated during exploratory analysis, so some expected filenames are still described within the individual scripts. The long-term goal is to route all inputs through the shared configuration file.
