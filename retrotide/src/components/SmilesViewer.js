import React from 'react';
import SmilesDrawerWrapper from './SmilesDrawerWrapper';

// This tab contains the text entry for writing SMILES format and includes
// a visualizer to display the molecule as entered, as well as check
// that the input string is a valid SMILES

class SmilesViewer extends React.Component {

  render() {
    return (
      <div className='Smiles form'>
        <p>Enter the molecule in SMILES format here</p>
        <SmilesDrawerWrapper />
      </div>
    )
  }
}


export default SmilesViewer;