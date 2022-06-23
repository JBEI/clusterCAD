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
      });
    }
    this.setState({DomainList: updatedDomainList}); 
  };

  // when a button is clicked, disable the other buttons
  // to avoid logical overlap errors
  toggleButtons = ClickedButtonName => {
    let updatedButtonList = this.state.ButtonList;

    // clicked button is active, ignore it
    // other buttons are toggled from previous
    // if we need some buttons to not toggle, we can add another
    // flag to them in future and filter those out here

    for (var ButtonKey in updatedButtonList) {
      if (updatedButtonList[ButtonKey].domainName !== ClickedButtonName) {
        updatedButtonList[ButtonKey].disabled = !updatedButtonList[ButtonKey].disabled;
      }
    }

    this.setState({ButtonList: updatedButtonList});
  }

  showOptionsModal = ClickedDomain => {
    if (ClickedDomain.options) {
      console.log(ClickedDomain.options);
    }
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
            {this.getAllButtons().map((DomainButton, index) => (
              <Button 
                className='addDomainButton' 
                disabled={DomainButton.disabled} 
                key={index} 
                onClick={ () => {
                  this.toggleDomains(DomainButton.domains); 
                  this.toggleButtons(DomainButton.domainName);
                } }
              >
                {DomainButton.domainName}
              </Button>
              ))
            }
          </div>
          <div className="DomainSandbox">
            {this.getPresentDomains().map((DomainDiv, index) => (
                <div key={DomainDiv.domainName + index} className="DomainWrapper" onClick={() => this.showOptionsModal(DomainDiv)}>
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