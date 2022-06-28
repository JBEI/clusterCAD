import React from 'react';
import { connect } from 'react-redux';
import { updateModuleDomains } from '../redux/actions/actions';
import Button from '../components/Button';

const mapDispatchToProps = dispatch => {
  return {
    updateModuleDomains: module => dispatch(updateModuleDomains(module)),
    dispatch,
  }
};

class ModuleBuilder extends React.Component {

  constructor(props) {
    super(props);
    this.state = {
      DomainList: props.domainList,
      ButtonList: props.buttonList,
      ModuleType: props.type,
      deleteFunction: props.deleteFunction,
      updateFunction: props.updateFunction,
      optionsModalOpen: false,
      optionsModalContent: {},
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
    this.props.updateFunction(this.props.id, updatedDomainList);
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
      this.setState({optionsModalContent: {domain: ClickedDomain.domainName, name: optionName, list: optionList, selected: selectedOption}});
      this.setState({optionsModalOpen: !this.state.optionsModalOpen});
    }
  }

  selectNewOption = (modalInfo, option) => {
    // this isn't really a list, it's an object
    let domains = this.state.DomainList;
    let {domain, name, list} = modalInfo;
    let updatedModule = {
      DomainList: {
        ...domains,
        [domain]: {
          ...domains[domain],
          options: {
            selected: option
      }}}
    }
    this.state.updateFunction(this.props.id, updatedModule);
    this.setState({DomainList: updatedModule});
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
            {this.state.optionsModalContent.name}
            :
            {this.state.optionsModalContent.list.map((option, index) => (
                <Button 
                  className={(option === this.state.optionsModalContent.selected) ? 'optionButton selected' : 'optionButton'}
                  key={'option'+index}
                  onClick={() => this.selectNewOption(this.state.optionsModalContent, option)}
                  >
                  {option}
                </Button>
              ))
            }
          </div>
        : null
        }        
      </div>
    )
  }

}

export default connect(null, mapDispatchToProps)(ModuleBuilder);
