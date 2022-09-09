import React from 'react';
import Button from '../components/Button';
import deleteIcon from '../images/delete-bin-7-fill.png';
import editIcon from '../images/edit-line.png';
import addIcon from '../images/add-circle-fill.png';

class ModuleBuilder extends React.Component {

  constructor(props) {
    super(props);
    this.state = {
      DomainList: props.DomainList,
      ButtonList: props.ButtonList,
      ModuleType: props.type,
      deleteFunction: props.deleteFunction,
      updateFunction: props.updateFunction,
      optionsModalOpen: false,
      optionsModalContent: {},
      buttonsModalOpen: false,
      isRemoving: false,
      isInserting: (props.type == 'extending' ? true : false),
    }
  };

  // this triggers the inserting module css animation
  // and adds a one second delay to allow the animation to run
  // note the iife is necessary to make setTimeout work
  componentDidMount() {
    setTimeout(() => this.setState({isInserting: false}), 1000);
  }

  // this populates the button list to add optional domains
  // it should maybe live in componentDidMount
  getAllButtons = () => {
    let allButtons = [];
    for (var ButtonObject in this.state.ButtonList) {
      allButtons.push(this.state.ButtonList[ButtonObject]);
    }
    return allButtons;
  };

  // this populates the domain sandbox with currently selected domains
  getPresentDomains = () => {
    let presentDomains = [];
    for (var DomainObject in this.state.DomainList) {
      if (this.state.DomainList[DomainObject].present) {
        presentDomains.push(this.state.DomainList[DomainObject]);
      } else {
        // we push a placeholder button that opens the buttons list to
        // insert optional modules
        presentDomains.push({button: "add"});
      }
    }
    return presentDomains;
  };

  // clicking the optional domain buttons adds or removes groups of domains
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

  // to insert domains, click the + in line and choose from the buttons in the modal
  // right now we only support one modal per module, so no props are needed
  // the buttons are already set up in another function so this is pretty basic
  toggleButtonsModal = () => {
    this.setState({buttonsModalOpen: !this.state.buttonsModalOpen});
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
      this.setState({optionsModalContent: {
                      domain: ClickedDomain.domainName, 
                      name: optionName, 
                      list: optionList, 
                      selected: selectedOption
                    }});
      this.setState({optionsModalOpen: !this.state.optionsModalOpen});
    }
  }

  // after opening the options modal, clicking an option button will update state
  selectNewOption = (modalInfo, option) => {
    // this isn't really a list, it's an object
    let currentOptions = this.state.optionsModalContent;
    let domains = this.state.DomainList;
    let {domain} = modalInfo;
    let updatedModule = {
      ...domains,
      [domain]: {
        ...domains[domain],
        options: {
          ...domains[domain].options,
          selected: option
      }}}
    this.state.updateFunction(this.props.id, updatedModule);
    this.setState({DomainList: updatedModule});
    // hacky doing this multiple times but it's not updating rn
    this.setState({optionsModalContent: {
      ...currentOptions,
      selected: option,
    }});
  }

  // this function triggers the css animation by setting the removing flag
  // then deleted the element after a short timeout
  // note the iife is necessary to make the timeout work
  deleteModule = () => {
    this.setState({isRemoving: true});
    setTimeout(() => this.state.deleteFunction(this.props.id), 1000);
  }

  // this takes the present domains and inserts either a domain (if required)
  // or a button that allows you to add optional domains 
  populateDomainSandbox = () => {
    return this.getPresentDomains().map((DomainDiv, index) => {
      if(DomainDiv.domainName) {
        console.log("found name");
        return (
          <div key={DomainDiv.domainName + index} className={`DomainWrapper ${DomainDiv.options ? 'MoreOptions' : ''}`} 
               onClick={() => this.toggleOptionsModal(DomainDiv)}>
            <div className={"Domain " + DomainDiv.domainName}>
              {DomainDiv.domainName}
            </div>
            <span className="optionsIcon">
              <img src={editIcon} />
            </span>
          </div>
        )
      } else {
        console.log("found button");
        return (
          <Button className="addDomains" onClick={() => {this.toggleButtonsModal()}}>
            <img src={addIcon} />
          </Button>
        )
      }
    });
  }

  render() {
    return (
      <div className={`ModuleBuilder 
        ${this.state.isInserting ? 'inserting' : ''} 
        ${this.state.isRemoving ? 'removing' : ''}`}
        >
        <div className="DomainHeader">
          <div className="DomainTitle">
            <div> Module {this.props.index + 1} </div>
            <div> {this.state.ModuleType} </div>
          </div>
          {this.state.ModuleType === 'extending' ? 
            <div className="DomainHeaderButton">
              <div className='DeleteModuleButton' onClick={() => {this.deleteModule()}}>
                <img src={deleteIcon} /> 
              </div> 
            </div>
            : null
          } 
        </div> 
        <div className="DomainToolbox">
          <div className={`DomainButtonList ${this.state.buttonsModalOpen ? 'open' : ''}`}>
            {this.getAllButtons().map((DomainButton, index) => (
              <Button 
                className='AddDomainButton' 
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
            { this.populateDomainSandbox() }
          </div>
        </div>
        {this.state.optionsModalOpen ? 
          <div className="OptionsModal">
            {this.state.optionsModalContent.name}
            :
            {this.state.optionsModalContent.list.map((option, index) => (
                <Button 
                  className={(option === this.state.optionsModalContent.selected) ? 'OptionButton selected' : 'OptionButton'}
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

export default ModuleBuilder;
