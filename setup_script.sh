#!/bin/bash

# Enhanced Setup script for the Protein Diversity Pipeline
# This script installs all required dependencies with improved error handling

set -e  # Exit on any error

echo "Protein Diversity Pipeline Setup"
echo "================================"
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print colored output
print_status() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_header() {
    echo -e "${BLUE}$1${NC}"
}

# Function to retry commands
retry_command() {
    local max_attempts=3
    local delay=5
    local attempt=1
    local cmd="$*"
    
    while [ $attempt -le $max_attempts ]; do
        print_status "Attempt $attempt/$max_attempts: $cmd"
        if eval "$cmd"; then
            return 0
        else
            if [ $attempt -lt $max_attempts ]; then
                print_warning "Command failed. Retrying in $delay seconds..."
                sleep $delay
                delay=$((delay * 2))  # Exponential backoff
            fi
            attempt=$((attempt + 1))
        fi
    done
    
    print_error "Command failed after $max_attempts attempts: $cmd"
    return 1
}

# Check if conda/mamba is available
check_conda() {
    if command -v mamba &> /dev/null; then
        CONDA_CMD="mamba"
        print_status "Found mamba (faster than conda)"
    elif command -v conda &> /dev/null; then
        CONDA_CMD="conda"
        print_status "Found conda"
    else
        print_error "Neither conda nor mamba found. Please install Miniconda/Anaconda first."
        echo "Download from: https://docs.conda.io/en/latest/miniconda.html"
        exit 1
    fi
}

# Create conda environment
# Create conda environment
create_environment() {
    print_header "Creating conda environment..."
    
    ENV_NAME="protein-pipeline"
    
    # Check if environment already exists
    if $CONDA_CMD env list | grep -q "^$ENV_NAME "; then
        print_warning "Environment '$ENV_NAME' already exists."
        read -p "Do you want to remove and recreate it? (y/N): " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            print_status "Removing existing environment..."
            $CONDA_CMD env remove -n $ENV_NAME -y
        else
            print_status "Using existing environment. Updating packages..."
            # Ensure environment is activated
            source $($CONDA_CMD info --base)/etc/profile.d/conda.sh
            $CONDA_CMD activate $ENV_NAME
            
            # Verify activation worked
            if [[ "$CONDA_DEFAULT_ENV" != "$ENV_NAME" ]]; then
                print_error "Failed to activate existing environment"
                exit 1
            fi
            return
        fi
    fi
    
    print_status "Creating new environment '$ENV_NAME'..."
    retry_command "$CONDA_CMD create -n $ENV_NAME python=3.9 -y"
    
    print_status "Activating environment..."
    source $($CONDA_CMD info --base)/etc/profile.d/conda.sh
    $CONDA_CMD activate $ENV_NAME
    
    # Verify activation worked
    if [[ "$CONDA_DEFAULT_ENV" != "$ENV_NAME" ]]; then
        print_error "Failed to activate conda environment '$ENV_NAME'"
        print_error "Please check your conda installation"
        exit 1
    fi
    
    print_status "✅ Environment '$ENV_NAME' created and activated"
    print_status "Current environment: $CONDA_DEFAULT_ENV"
}

# Verify conda environment is active
verify_environment() {
    if [[ -z "$CONDA_DEFAULT_ENV" ]]; then
        print_error "No conda environment is active"
        print_error "This should not happen. Please check the environment activation."
        exit 1
    fi
    
    if [[ "$CONDA_DEFAULT_ENV" != "protein-pipeline" ]]; then
        print_error "Wrong conda environment active: $CONDA_DEFAULT_ENV"
        print_error "Expected: protein-pipeline"
        exit 1
    fi
    
    print_status "✅ Conda environment verified: $CONDA_DEFAULT_ENV"
}

# Install Python packages with minimal GUI dependencies
install_python_packages() {
    print_header "Installing Python packages..."
    
    # Verify environment is active
    verify_environment
    
    # Configure conda for better performance and reliability
    print_status "Configuring conda settings..."
    $CONDA_CMD config --set channel_priority strict
    $CONDA_CMD config --set remote_read_timeout_secs 120
    
    # Add conda-forge channel (but skip bioconda to avoid conflicts)
    print_status "Adding conda-forge channel..."
    $CONDA_CMD config --add channels conda-forge 2>/dev/null || true
    
    # Install core packages without heavy GUI dependencies first
    print_status "Installing essential scientific packages..."
    retry_command "$CONDA_CMD install -c conda-forge numpy scipy pandas -y"
    
    print_status "Installing machine learning packages..."
    retry_command "$CONDA_CMD install -c conda-forge scikit-learn joblib threadpoolctl -y"
    
    print_status "Installing utilities..."
    retry_command "$CONDA_CMD install -c conda-forge tqdm -y"
    
    # Install matplotlib with minimal backend
    print_status "Installing matplotlib (minimal backend)..."
    retry_command "$CONDA_CMD install -c conda-forge matplotlib-base -y"
    
    # Install seaborn (depends on matplotlib)
    print_status "Installing seaborn..."
    retry_command "$CONDA_CMD install -c conda-forge seaborn-base -y"
    
    # Install pip packages (more reliable for biopython)
    print_status "Upgrading pip in conda environment..."
    retry_command "python -m pip install --upgrade pip"
    
    print_status "Installing BioPython via pip..."
    retry_command "python -m pip install biopython"
    
    print_status "Installing additional packages via pip..."
    retry_command "python -m pip install requests urllib3"
    
    # Set matplotlib to use non-GUI backend
    print_status "Configuring matplotlib backend..."
    python -c "
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
plt.ioff()  # Turn off interactive mode
print('✅ Matplotlib configured for non-interactive use')
"
    
    print_status "✅ All Python packages installed successfully"
}

# Install MUSCLE
install_muscle() {
    print_header "Installing MUSCLE for multiple sequence alignment..."
    
    # Try conda installation first
    print_status "Installing MUSCLE v5 via conda..."
    if retry_command "$CONDA_CMD install -c bioconda muscle -y"; then
        print_status "MUSCLE installed successfully via conda"
        
        # Test MUSCLE installation
        if muscle -version &>/dev/null; then
            print_status "✅ MUSCLE is working correctly"
            return 0
        else
            print_warning "MUSCLE installed but version check failed, trying binary installation..."
        fi
    else
        print_warning "Conda installation failed, trying binary installation..."
    fi
    
    # Alternative: download MUSCLE binary
    print_status "Installing MUSCLE binary..."
    case "$(uname -s)" in
        Linux*)
            MUSCLE_URL="https://github.com/rcedgar/muscle/releases/download/v5.1/muscle5.1.linux_intel64"
            ;;
        Darwin*)
            MUSCLE_URL="https://github.com/rcedgar/muscle/releases/download/v5.1/muscle5.1.macos_intel64"
            ;;
        *)
            print_error "Unsupported operating system for MUSCLE binary download"
            return 1
            ;;
    esac
    
    print_status "Downloading MUSCLE binary from GitHub..."
    if command -v wget &> /dev/null; then
        retry_command "wget -O muscle '$MUSCLE_URL'"
    elif command -v curl &> /dev/null; then
        retry_command "curl -L -o muscle '$MUSCLE_URL'"
    else
        print_error "Neither wget nor curl found. Cannot download MUSCLE binary."
        return 1
    fi
    
    chmod +x muscle
    
    # Move to a directory in PATH
    if [[ -w "/usr/local/bin" ]] && sudo -n true 2>/dev/null; then
        sudo mv muscle /usr/local/bin/
        print_status "MUSCLE installed to /usr/local/bin/"
    else
        mkdir -p ~/.local/bin
        mv muscle ~/.local/bin/
        export PATH="$HOME/.local/bin:$PATH"
        
        # Add to shell configuration
        for rc_file in ~/.bashrc ~/.zshrc; do
            if [[ -f "$rc_file" ]]; then
                if ! grep -q 'export PATH="$HOME/.local/bin:$PATH"' "$rc_file"; then
                    echo 'export PATH="$HOME/.local/bin:$PATH"' >> "$rc_file"
                fi
            fi
        done
        print_status "MUSCLE installed to ~/.local/bin/"
    fi
    
    # Test installation
    if command -v muscle &> /dev/null && muscle -version &>/dev/null; then
        print_status "✅ MUSCLE binary installed and working"
    else
        print_error "❌ MUSCLE installation failed"
        return 1
    fi
}

# Install DSSP
install_dssp() {
    print_header "Installing DSSP for Structure Analysis..."
    # Try conda installation first
    print_status "Installing DSSP via conda..."
    if retry_command "$CONDA_CMD install -c conda-forge dssp -y"; then
        print_status "✅ DSSP installed successfully"
    else
        print_error "❌ DSSP installation failed via conda"
        print_status "You can try installing it manually with: sudo apt install dssp"
        return 1
    fi
}

# Test installation
test_installation() {
    print_header "Testing installation..."
    
    # Test Python packages
    print_status "Testing Python packages..."
    python -c "
import numpy, scipy, sklearn, Bio
from Bio import SeqIO, AlignIO
from Bio.Align import PairwiseAligner
print('✅ All Python packages imported successfully')
"
    
    # Test MUSCLE
    print_status "Testing MUSCLE..."
    if command -v muscle &> /dev/null; then
        muscle -version
        print_status "✅ MUSCLE is working"
    else
        print_warning "❌ MUSCLE not found in PATH"
    fi
    
    # Test DSSP
    print_status "Testing DSSP..."
    if command -v mkdssp &> /dev/null; then
        mkdssp -h > /dev/null 2>&1
        print_status "✅ DSSP is working"
    else
        print_warning "❌ DSSP not found in PATH"
    fi
}


# Main installation function
main() {
    print_header "Protein Diversity Pipeline Setup"
    echo ""
    
    # Check system requirements
    check_conda
    
    # Install everything
    create_environment
    install_python_packages
    install_muscle
    install_dssp
    
    print_header "✅ Setup completed successfully!"
    echo ""
    print_status "NEXT STEPS:"
    echo "1. Activate the environment: conda activate protein-pipeline"
    echo "2. Test the installation: python ProteinDiversity_pipeline.py --check-requirements"
    echo "3. Run with test data: python ProteinDiversity_pipeline.py test_proteins.fasta"
    echo "4. Use your own data: python ProteinDiversity_pipeline.py your_proteins.fasta"
    echo ""
    print_status "IMPORTANT NOTES:"
    echo "• Always activate the environment before running the pipeline"
    echo "• Check the comprehensive reports in the output directory"
    echo ""
    print_status "TROUBLESHOOTING:"
    echo "• If MUSCLE fails: conda install -c bioconda muscle=5"
    echo "• If DSSP fails: sudo apt install dssp"
    echo "• For Python package issues: Ensure you are in the correct conda environment"
    echo "• For BioPython issues: Check if Biopython is installed in the environment"
    echo "• For matplotlib issues: Ensure you are using the non-interactive backend"
    echo "• For seaborn issues: Ensure seaborn is installed in the environment"
    echo "• For conda/mamba issues: Ensure conda/mamba is installed and available in PATH"
    echo "• For environment activation issues: Check if activate_pipeline.sh is sourced correctly"
    echo "• For general issues: Check the logs in the output directory for more details"
}

# Handle command line arguments
case "${1:-}" in
    --help|-h)
        echo "Protein Diversity Pipeline Setup Script"
        echo ""
        echo "Usage: $0 [OPTIONS]"
        echo ""
        echo "Options:"
        echo "  --help, -h          Show this help message"
        echo "  --test-only         Only run tests, skip installation"
        echo "  --create-env-only   Only create environment, skip other steps"
        echo ""
        echo "This script will:"
        echo "1. Create a conda environment with all dependencies"
        echo "2. Install MUSCLE for multiple sequence alignment"
        echo "3. Install DSSP for structure analysis"
        exit 0
        ;;
    --test-only)
        test_installation
        exit 0
        ;;
    --create-env-only)
        check_conda
        create_environment
        install_python_packages
        exit 0
        ;;
    *)
        main
        ;;
esac

