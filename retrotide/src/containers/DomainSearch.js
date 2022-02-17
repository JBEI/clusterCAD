import React from 'react';
import Button from '../components/Button';
import Draggable from 'react-draggable';

class DomainSearch extends React.Component {

  constructor(props) {
    super(props);
    this.state = {modules: []}
  }

  insertDomain = d => {
    this.setState({modules: [...this.state.modules, d]});
    console.log(this.state.modules);
  }

  deleteDomain = d => {
    let dindex = d.index;
    let oldModules = this.state.modules;
    this.setState({modules: oldModules.slice(0, dindex).concat(dindex+1, oldModules.length())});
  }

  render() {
    return (
      <div className='DomainSearch'>
        <h3>Construct Module</h3>
        <div className="DomainToolbox">
          <div className="DomainButtonList">
            <Button onClick={ () => {this.insertDomain("A")} }> A
            </Button>
          </div>
          <div className="DomainSandbox">
            {this.state.modules.map(M => (
              <Draggable>
                <div className={"Domain handle " + M.value}>
                A
                </div>
              </Draggable>
              ))
            }
          </div>
        </div>
      </div>
    )
  }
};

export default DomainSearch;