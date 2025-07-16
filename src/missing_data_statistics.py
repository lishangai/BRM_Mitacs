
"""
Missing Data Statistics Script
Statistics of zero values proportion in each row and column of CSV files
"""

import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Define Project Root for robust path management
PROJECT_ROOT = Path(__file__).resolve().parent.parent

# Set font for Chinese characters (if needed)
plt.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei']
plt.rcParams['axes.unicode_minus'] = False

def calculate_zero_proportion(data, row_names=None):
    """
    Calculate proportion of zero values in data
    
    Parameters:
    -----------
    data : pandas.DataFrame
        Input data
    row_names : pandas.Series, optional
        Row names/identifiers
        
    Returns:
    --------
    row_zero_prop : pandas.Series
        Zero proportion for each row
    col_zero_prop : pandas.Series  
        Zero proportion for each column
    """
    # Calculate zero proportion for each row
    row_zero_counts = (data == 0).sum(axis=1)
    row_zero_prop = row_zero_counts / data.shape[1]
    
    # If row names are provided, use them as index
    if row_names is not None:
        row_zero_prop.index = row_names
    
    # Calculate zero proportion for each column
    col_zero_counts = (data == 0).sum(axis=0)
    col_zero_prop = col_zero_counts / data.shape[0]
    
    return row_zero_prop, col_zero_prop

def plot_zero_proportions(row_zero_prop, col_zero_prop, output_dir="./"):
    """
    Plot zero proportion distributions
    
    Parameters:
    -----------
    row_zero_prop : pandas.Series
        Zero proportion for each row
    col_zero_prop : pandas.Series
        Zero proportion for each column
    output_dir : str
        Output directory
    """
    # Create output directory
    Path(output_dir).mkdir(exist_ok=True)
    
    # Plot row zero proportion distribution
    plt.figure(figsize=(12, 5))
    
    plt.subplot(1, 2, 1)
    plt.hist(row_zero_prop, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
    plt.xlabel('Zero Proportion per Row')
    plt.ylabel('Frequency')
    plt.title('Distribution of Zero Proportions per Row')
    plt.grid(True, alpha=0.3)
    
    # Plot column zero proportion distribution
    plt.subplot(1, 2, 2)
    plt.hist(col_zero_prop, bins=50, alpha=0.7, color='lightcoral', edgecolor='black')
    plt.xlabel('Zero Proportion per Column')
    plt.ylabel('Frequency')
    plt.title('Distribution of Zero Proportions per Column')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/zero_proportion_distribution.png", dpi=300, bbox_inches='tight')
    plt.show()
    
    # Plot trend lines (if data is not too large)
    if len(col_zero_prop) <= 500 and len(row_zero_prop) <= 500:
        plt.figure(figsize=(10, 8))
        
        plt.subplot(2, 1, 1)
        plt.plot(range(len(col_zero_prop)), col_zero_prop, 'b-', alpha=0.7, linewidth=1)
        plt.xlabel('Column Index')
        plt.ylabel('Zero Proportion')
        plt.title('Zero Proportion Trend across Columns')
        plt.grid(True, alpha=0.3)
        
        plt.subplot(2, 1, 2)
        plt.plot(range(len(row_zero_prop)), row_zero_prop, 'r-', alpha=0.7, linewidth=1)
        plt.xlabel('Row Index')
        plt.ylabel('Zero Proportion')
        plt.title('Zero Proportion Trend across Rows')
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/zero_proportion_trends.png", dpi=300, bbox_inches='tight')
        plt.show()

def generate_summary_report(data, row_zero_prop, col_zero_prop, output_dir="./", row_names=None):
    """
    Generate summary report
    
    Parameters:
    -----------
    data : pandas.DataFrame
        Original data
    row_zero_prop : pandas.Series
        Zero proportion for each row
    col_zero_prop : pandas.Series
        Zero proportion for each column
    output_dir : str
        Output directory
    row_names : pandas.Series, optional
        Row names/identifiers
    """
    # Create output directory
    Path(output_dir).mkdir(exist_ok=True)
    
    # Generate summary statistics
    summary = {
        "Basic Data Information": {
            "Total Rows": data.shape[0],
            "Total Columns": data.shape[1],
            "Total Elements": data.size,
            "Total Zero Values": (data == 0).sum().sum(),
            "Overall Zero Proportion": (data == 0).sum().sum() / data.size
        },
        "Row Zero Proportion Statistics": {
            "Minimum": row_zero_prop.min(),
            "Maximum": row_zero_prop.max(),
            "Mean": row_zero_prop.mean(),
            "Median": row_zero_prop.median(),
            "Standard Deviation": row_zero_prop.std(),
            "Completely Zero Rows": (row_zero_prop == 1.0).sum(),
            "Completely Non-zero Rows": (row_zero_prop == 0.0).sum()
        },
        "Column Zero Proportion Statistics": {
            "Minimum": col_zero_prop.min(),
            "Maximum": col_zero_prop.max(),
            "Mean": col_zero_prop.mean(),
            "Median": col_zero_prop.median(),
            "Standard Deviation": col_zero_prop.std(),
            "Completely Zero Columns": (col_zero_prop == 1.0).sum(),
            "Completely Non-zero Columns": (col_zero_prop == 0.0).sum()
        }
    }
    
    # Print summary report
    print("\n" + "="*60)
    print("           Missing Data (Zero Values) Statistics Report")
    print("="*60)
    
    for category, stats in summary.items():
        print(f"\n{category}:")
        print("-" * 40)
        for key, value in stats.items():
            if isinstance(value, float):
                print(f"  {key}: {value:.4f}")
            else:
                print(f"  {key}: {value}")
    
    # Save detailed results to CSV
    # For rows - use row names if available
    if row_names is not None:
        row_results = pd.DataFrame({
            'row_name': row_names.values if hasattr(row_names, 'values') else row_names,
            'zero_proportion': row_zero_prop.values,
            'zero_count': (data == 0).sum(axis=1).values,
            'total_count': data.shape[1]
        })
    else:
        row_results = pd.DataFrame({
            'row_index': range(len(row_zero_prop)),
            'zero_proportion': row_zero_prop.values,
            'zero_count': (data == 0).sum(axis=1).values,
            'total_count': data.shape[1]
        })
    row_results.to_csv(f"{output_dir}/row_zero_statistics.csv", index=False, encoding='utf-8-sig')
    
    # For columns
    col_results = pd.DataFrame({
        'column_name': data.columns if hasattr(data, 'columns') else range(len(col_zero_prop)),
        'zero_proportion': col_zero_prop.values,
        'zero_count': (data == 0).sum(axis=0).values,
        'total_count': data.shape[0]
    })
    col_results.to_csv(f"{output_dir}/column_zero_statistics.csv", index=False, encoding='utf-8-sig')
    
    # Save summary report
    with open(f"{output_dir}/summary_report.txt", 'w', encoding='utf-8') as f:
        f.write("Missing Data (Zero Values) Statistics Report\n")
        f.write("="*60 + "\n\n")
        
        for category, stats in summary.items():
            f.write(f"{category}:\n")
            f.write("-" * 40 + "\n")
            for key, value in stats.items():
                if isinstance(value, float):
                    f.write(f"  {key}: {value:.4f}\n")
                else:
                    f.write(f"  {key}: {value}\n")
            f.write("\n")
    
    print(f"\nDetailed results saved to:")
    print(f"  - Row statistics: {output_dir}/row_zero_statistics.csv")
    print(f"  - Column statistics: {output_dir}/column_zero_statistics.csv") 
    print(f"  - Summary report: {output_dir}/summary_report.txt")

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Calculate zero value proportions for each row and column in CSV files')
    parser.add_argument('input_file', help='Input CSV file path (relative to project root)')
    parser.add_argument('--output_dir', '-o', default=None, help='Output directory. If not specified, a folder with the script name inside results/ will be used.')
    parser.add_argument('--delimiter', '-d', default=',', help='CSV delimiter (default: comma)')
    parser.add_argument('--header', action='store_true', help='Whether CSV file contains header')
    parser.add_argument('--encoding', '-e', default='utf-8', help='File encoding (default: utf-8)')
    parser.add_argument('--sample_size', '-s', type=int, help='Sample size (for large files)')
    parser.add_argument('--row_names_col', '-r', type=int, default=0, help='Column index for row names (default: 0, first column)')
    
    args = parser.parse_args()
    
    # Determine output directory based on convention, relative to project root
    if args.output_dir:
        # If user provides a path, resolve it relative to the project root
        output_dir = PROJECT_ROOT / args.output_dir
    else:
        # Default to <PROJECT_ROOT>/results/<script_name>/
        script_name = Path(__file__).stem
        output_dir = PROJECT_ROOT / "results" / script_name
        
    try:
        # Construct absolute path for the input file from project root
        input_file_path = PROJECT_ROOT / args.input_file
        print(f"Reading file: {input_file_path}")
        
        # Read CSV file
        read_kwargs = {
            'delimiter': args.delimiter,
            'encoding': args.encoding,
            'header': 0 if args.header else None
        }
        
        if args.sample_size:
            # Sample large files
            data = pd.read_csv(input_file_path, nrows=args.sample_size, **read_kwargs)
            print(f"Sampled first {args.sample_size} rows")
        else:
            data = pd.read_csv(input_file_path, **read_kwargs)
        
        print(f"Data shape: {data.shape}")
        
        # Extract row names if specified
        row_names = None
        if args.row_names_col is not None and args.row_names_col < data.shape[1]:
            row_names = data.iloc[:, args.row_names_col]
            # Remove the row names column from data for analysis
            data_for_analysis = data.drop(data.columns[args.row_names_col], axis=1)
            print(f"Using column {args.row_names_col} ({data.columns[args.row_names_col]}) as row names")
        else:
            data_for_analysis = data
        
        # Ensure data is numeric
        numeric_columns = data_for_analysis.select_dtypes(include=[np.number]).columns
        if len(numeric_columns) < data_for_analysis.shape[1]:
            print(f"Note: Detected non-numeric columns, will only analyze numeric columns ({len(numeric_columns)}/{data_for_analysis.shape[1]})")
            data_for_analysis = data_for_analysis[numeric_columns]
        
        # Calculate zero proportions
        print("Calculating zero proportions...")
        row_zero_prop, col_zero_prop = calculate_zero_proportion(data_for_analysis, row_names)
        
        # Generate summary report
        generate_summary_report(data_for_analysis, row_zero_prop, col_zero_prop, output_dir, row_names)
        
        # Plot charts
        print("Generating visualization charts...")
        plot_zero_proportions(row_zero_prop, col_zero_prop, output_dir)
        
        print("\nAnalysis completed!")
        
    except Exception as e:
        print(f"Error: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    exit(main()) 