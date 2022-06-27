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
      optionsModalOpen: false,
      optionsModalContent: [],
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
            ...selectedDomain,
            present: false,
          }
          updatedDomainList = {
            ...updatedDomainList,
            [Domain]: deleteDomain,
          }
        } else {
          let insertDomain = {
            ...selectedDomain,
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

  // some domains (AT and KR) have settable properties
  // we have 2 tuples, the currently selected option, and the name and list of options
  toggleOptionsModal = ClickedDomain => {
    if (ClickedDomain.options) {
      let optionName;
      let optionList;
      let selectedOption;
      for (var option in ClickedDomain.options) {
        if (option === 'selected') {
          selectedOption = ClickedDomain.options[option];
        } else {
          optionName = option;
          optionList = ClickedDomain.options[option];
        }
      }
      this.setState({optionsModalContent: optionName});
      this.setState({optionsModalOpen: !this.state.optionsModalOpen});
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
                <div key={DomainDiv.domainName + index} className="DomainWrapper" onClick={() => this.toggleOptionsModal(DomainDiv)}>
                  <div className={"Domain " + DomainDiv.domainName}>
                    {DomainDiv.domainName}
                  </div>
                </div>
              ))
            }
          </div>
        </div>
        {this.state.optionsModalOpen ? 
          <div className="optionsModal">
            {this.state.optionsModalContent}
            :
            {/* buttons with one selected and onclick to select */}
            A B C D
          </div>
        : null
        }        
      </div>
    )
  }

}

export default ModuleBuilder;