import React from 'react';
import Button from '../components/Button';

class ModuleBuilder extends React.Component {

  // domain list gets passed
  constructor(props) {
    super(props);
    this.state = {
      DomainList: props.domainList,
      ModuleType: 'loading',
    }
  }

  getAllDomains = () => {
    let allDomains = [];
    for (var DomainObject in this.state.DomainList) {
      allDomains.push(this.state.DomainList[DomainObject]);
    }
    return allDomains;
  };

  getPresentDomains = () => {
    let presentDomains = [];
    for (var DomainObject in this.state.DomainList) {
      if (this.state.DomainList[DomainObject].present) {
        presentDomains.push(this.state.DomainList[DomainObject]);
      }
    }
    return presentDomains;
  };

  insertDomain = Domain => {
    let selectedDomain = this.state.DomainList[Domain];

    if(selectedDomain.present) {
      console.log("ERR that domain is already selected " + selectedDomain.domainName);
      return -1;
    } else {
      // we'll need some logic here to add multiple nodes depending on selected name
      let insertDomain = {
        domainName: Domain,
        present: true,
      }
      let updatedDomainList = {
        ...this.state.DomainList,
        [Domain]: insertDomain,
      };
      this.setState({DomainList: updatedDomainList});      
    }
  }

  deleteDomain = Domain => {
    let deleteDomain = {
      domainName: Domain,
      present: false,
    }
    let updatedDomainList = {
      ...this.state.DomainList,
      [Domain]: deleteDomain,
    };
    this.setState({DomainList: updatedDomainList});
  }

  render() {
    return (
      <div className='ModuleBuilder'>
        <div className="DomainToolbox">
          <div className="DomainButtonList">
            <select className="ModuleType">
              <option value="loading">Loading</option>
              <option value="extending">Extending</option>
              <option value="terminating">Terminating</option>
            </select>
            {this.getAllDomains().map((DomainButton, index) => (
              <Button className='addDomainButton' key={index} onClick={ () => {this.insertDomain(DomainButton.domainName)} }>
                {DomainButton.domainName}
              </Button>
              ))
            }
          </div>
          <div className="DomainSandbox">
            {this.getPresentDomains().map((DomainDiv, index) => (
                <div key={DomainDiv.domainName + index}>
                  <div className={"Domain " + DomainDiv.domainName}>
                    {DomainDiv.domainName}
                    <div className="deleteIcon" onClick={ () => {this.deleteDomain(DomainDiv.domainName)} }> X </div>
                  </div>
                </div>
              ))
            }
          </div>
        </div>
      </div>
    )
  }

}

export default ModuleBuilder;