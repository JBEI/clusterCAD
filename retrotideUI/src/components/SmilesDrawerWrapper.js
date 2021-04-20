import React from 'react';
import Button from './Button';
import { connect } from 'react-redux';
import SmilesDrawer from 'smiles-drawer';
import { updateMolecule } from '../redux/actions/actions';

// this wrapper wraps together an input and connected canvas into one unit
// in future we may want to break out the input to make things more versatile

let canvasOptions = {
  height: 500,
  width: 800,
};

const mapDispatchToProps = dispatch => {
  return {
    updateMolecule: molecule => dispatch(updateMolecule(molecule)),
    dispatch,
  }
};

class SmilesDrawerWrapper extends React.Component {

  constructor(props) {
    super(props);
    this.state = {
      errorMessage: "",
      smilesString: ""
    }
    this.inputRef = React.createRef(); // drawer needs access to the input node
    this.errorRef = React.createRef(); // notifications of invaid input appear here
    this.canvasRef = React.createRef(); // drawer needs to know where to draw
  }

  // put local state into global state
  submitMolecule = () => {
    let molecule = this.state.smilesString;
    // dispatch update action with local state value
    this.props.updateMolecule(molecule);
  }

  resetCanvas = () => {
    this.setState({
      smilesString: "",
      errorMessage: ""
    });
    this.canvasRef.current.getContext('2d').clearRect(0,0,2000,2000);
    this.inputRef.current.value = "";
  }

  parseStringInput = () => {
    let mol = this.inputRef.current.value; // we're double-parsing due to workaround
    SmilesDrawer.parse(mol,
      (tree) => {
        console.log(tree);
        this.setState({
          errorMessage: "",
          smilesString: mol
        }, () => {
          SmilesDrawer.apply(canvasOptions); // callback compensates for async setState
        });
      },
      (err) => {
        console.log(err);
        this.setState({
          errorMessage: err.name + " : " + err.message
        });
        console.log(this.state.errorMessage);
      }
    );
  }

  render() {
    let disableSubmit = (this.state.errorMessage.length > 0 || this.state.smilesString.length < 1);

    return (
      <div className='toolWrapper'>
        <input className="smilesInput" ref={this.inputRef} onChange={this.parseStringInput} type="text" />
        <Button onClick={this.submitMolecule} disabled={disableSubmit} className='submit'>Submit Molecule</Button>
        <Button onClick={this.resetCanvas} className='clear'>Clear Canvas</Button>
        <p className='errorField' ref={this.errorRef} >{this.state.errorMessage}</p> {/* any parsing errors will appear here */}
        <canvas className="smilesCanvas" width={canvasOptions.width} height={canvasOptions.height} ref={this.canvasRef} data-smiles={this.state.smilesString}></canvas>
        <p className="smilesCanvasAcknowledgment"> Smiles Drawer is used courtesy of Daniel Probst and Jean-Louis Reymond 
          (<a href="http://dx.doi.org/10.1021/acs.jcim.7b00425" target="_blank" rel="noopener noreferrer">original paper</a>)
        </p>
      </div>
    )
  }
}

export default connect(null, mapDispatchToProps)(SmilesDrawerWrapper);