import React from 'react';
import Button from '../components/Button';

class ModuleBuilder extends React.Component {

  // domain list gets passed
  constructor(props) {
    super(props);
    this.state = {
      DomainList: props.domainList,
      ModuleType: props.type,
      deleteFunction: props.deleteFunction,
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

  insertDomains = NewDomains => {
    for(var domain in NewDomains) {
      this.insertDomain(domain);
    }
  }

  toggleDomain = Domain => {
    let selectedDomain = this.state.DomainList[Domain];
    let updatedDomainList;

    if(selectedDomain.present) {
      let deleteDomain = {
        domainName: Domain,
        present: false,
      }
      updatedDomainList = {
        ...this.state.DomainList,
        [Domain]: deleteDomain,
      };
    } else {
      // we'll need some logic here to add multiple nodes depending on selected name
      let insertDomain = {
        domainName: Domain,
        present: true,
      }
      updatedDomainList = {
        ...this.state.DomainList,
        [Domain]: insertDomain,
      };     
    }
    this.setState({DomainList: updatedDomainList}); 
  }

  deleteDomain = Domain => {
    let domainToDelete = {
      domainName: Domain,
      present: false,
    }
    let updatedDomainList = {
      ...this.state.DomainList,
      [Domain]: domainToDelete,
    };
    this.setState({DomainList: updatedDomainList});
  }

  render() {
    return (
      <div className='ModuleBuilder'>
        <div className="DomainHeader">
          <div> Module {this.props.index + 1} </div>
          <div> {this.state.ModuleType} </div>
        </div>
        {this.state.ModuleType === 'extending' ? 
          <div className="DomainHeaderButton">
            <Button className='deleteModuleButton' onClick={() => {this.state.deleteFunction(this.props.id)}}> X </Button> 
          </div>
          : null
        }        
        <div className="DomainToolbox">
          <div className="DomainButtonList">
            {this.getAllDomains().map((DomainButton, index) => (
              <Button className='addDomainButton' key={index} onClick={ () => {this.toggleDomain(DomainButton.domainName)} }>
                {DomainButton.domainName}
              </Button>
              ))
            }
          </div>
          <div className="DomainSandbox">
            {this.getPresentDomains().map((DomainDiv, index) => (
                <div key={DomainDiv.domainName + index} className="DomainWrapper" >
                  <div className={"Domain " + DomainDiv.domainName}>
                    {DomainDiv.domainName}
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