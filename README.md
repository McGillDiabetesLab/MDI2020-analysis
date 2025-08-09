
# MDI2020-analysis

This repository contains MATLAB code for analyzing data from the MDI2020 Study, organized into several main components:

## Repository Structure

- **@MAPData/**: Core data handling and analysis (glucose, insulin, meals, outcomes).
- **@MAPPatient/**: Patient-level data and metadata management.
- **@MAPStudy/**: Study-level management, including multiple patients and interventions.
- **@MAPUtils/**: Utility functions for data manipulation, plotting, and statistics.
- **Data/**: Contains sample data files (e.g., `MDI2020.mat`).
- **exportData.m**: Main entry.

**Typical Workflow:**
1. Place your data in the `Data` folder.
2. Use the classes in `@MAPData`, `@MAPPatient`, and `@MAPStudy` to load, process, and analyze study data.
3. Analyze glucose, insulin, and meal outcomes with built-in methods.
4. Generate figures and tables using provided methods (e.g., `toFigure`, `toSheet`).
5. Export results to Excel, CSV, or MAT files for reporting.
6. Visualize data and summary statistics using plotting functions.

## Technical Overview: Data Export Workflow

### Main Entry: `exportData.m`
The main entry point for exporting study results is `exportData.m`. This script:
- Loads the study data from `Data/MDI2020.mat` using the `MAPStudy` class.
- Renames data fields for clarity.
- Calls `sMAP.toExcel('summary.xlsx', ...)` to export results to an Excel file.
- Runtime of less than 1 minute
- Generate a file: summary.xlsx

### Export Logic: `toExcel` and `toSheet`
The export process is handled by the `toExcel` method in `MAPUtils`, which:
1. Ensures the output folder and filename are set.
2. Calls `toSheet` to generate a structure of tables (sheets) summarizing study, patient, and outcome data.
3. Iterates over each sheet and writes it to the Excel file using MATLAB's `writetable` function.
4. Customizes the Excel file formatting for readability.

#### Pseudocode for `toExcel`:
```
function toExcel(obj, filename, ...)
	Ensure output folder and filename
	sheets = obj.toSheet(...)
	for each sheet in sheets:
		Write sheet to Excel
	Customize Excel formatting
end
```

#### How `toSheet` Works:
- Aggregates results for demographics, outcomes, comparisons, and individual patients.
- Supports options for study design (matched, crossover), summary statistics, and day/night intervals.
- Returns a struct where each field is a table to be exported as a separate Excel sheet.

This modular workflow allows flexible export of all relevant study results for reporting and further analysis.

## Contact

For questions or contributions, please open an issue or contact the repository maintainers.
