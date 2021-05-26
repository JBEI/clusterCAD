import React from 'react';
import MoleculeDiagram from '../components/MoleculeDiagram';
import SmilesViewer from '../components/SmilesViewer';
import Button from '../components/Button';

// all of this needs to move to presentational later but we'll see
// if we end up using that paradigm
class Design extends React.Component {

  constructor(props) {
    super(props);
    this.state = {entryMethod: 'smiles'}
    console.log(this.state);
  }

  render() {
    return (
      <div className='DesignContents'>
        <h3>Design PKS</h3>
        <div className="inputSelector">
          <Button className={"smilesButton" + (this.state.entryMethod === 'smiles' ? " selected" : "")}
            onClick={() => this.setState({entryMethod: 'smiles'})}
            >SMILES
          </Button>
          <Button 
            className={"drawerButton" + (this.state.entryMethod === 'drawMolecule' ? " selected" : "")}
            onClick={() => this.setState({entryMethod: 'drawMolecule'})}
            >Draw Molecule
          </Button>
        </div>
        <div>
          {this.state.entryMethod === 'smiles' ? 
            <SmilesViewer /> :
            <MoleculeDiagram />
          }
        </div>
      </div>
    )
  }
};

export default Design;