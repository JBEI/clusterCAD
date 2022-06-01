import React from 'react';
import Button from '../components/Button';

class ModuleBuilder extends React.Component {

  constructor(props) {
    super(props);
    this.state = {
      DomainList: props.domainList,
      ButtonList: props.buttonList,
      ModuleType: props.type,
      deleteFunction: props.deleteFunction,
    }
    console.log(this.state.ButtonList);
  };

  getAllButtons = () => {
    let allButtons = [];
    for (var ButtonObject in this.state.ButtonList) {
      allButtons.push(this.state.ButtonList[ButtonObject]);
    }
    return allButtons;
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
  };

  toggleDomains = domainsToToggle => {
    let updatedDomainList = this.state.DomainList;

    if (domainsToToggle.length > 0) {
      domainsToToggle.forEach((Domain) => {
        let selectedDomain = this.state.DomainList[Domain];

        if(selectedDomain.present) {
          let deleteDomain = {
            domainName: Domain,
            present: false,
          }
          updatedDomainList = {
            ...updatedDomainList,
            [Domain]: deleteDomain,
          }
        } else {
          let insertDomain = {
            domainName: Domain,
            present: true,
          }
          updatedDomainList = {
            ...updatedDomainList,
            [Domain]: insertDomain,
          }   
        }
        console.log(updatedDomainList);
      });
    }
    console.log(updatedDomainList);
    this.setState({DomainList: updatedDomainList}); 
  };

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
            {this.getAllButtons().map((DomainButton, index) => (
              <Button className='addDomainButton' key={index} onClick={ () => {this.toggleDomains(DomainButton.domains)} }>
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