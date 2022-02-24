import React from 'react';
import Button from '../components/Button';
import Draggable from 'react-draggable';

class DomainSearch extends React.Component {

  constructor(props) {
    super(props);
    this.state = {
      DomainList: {
        ACP: {domainName: 'ACP', present: false, index: -1}, 
        AT: {domainName: 'AT', present: false, index: -1},
        DH: {domainName: 'DH', present: false, index: -1},
        ER: {domainName: 'ER', present: false, index: -1},
        KR: {domainName: 'KR', present: false, index: -1},
        KS: {domainName: 'KS', present: false, index: -1},
      },
      insertIndex: 0,
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

  // getPresentDomainsInOrder

  insertDomain = Domain => {
    let selectedDomain = this.state.DomainList[Domain];
    let currentIndex = this.state.insertIndex;

    if(selectedDomain.present) {
      console.log("ERR that domain is already selected " + selectedDomain.domainName);
      return -1;
    } else {
      // we'll need some logic here to add multiple nodes depending on selected name
      let insertDomain = {
        domainName: Domain,
        present: true,
        index: currentIndex+1,
      }
      let updatedDomainList = {
        ...this.state.DomainList,
        [Domain]: insertDomain,
      };
      this.setState({insertIndex: currentIndex++});
      this.setState({DomainList: updatedDomainList});
      return currentIndex+1; // fix this
    }
  }

  deleteDomain = Domain => {
    let currentIndex = this.state.insertIndex;
      let deleteDomain = {
        domainName: Domain,
        present: false,
        index: -1,
      }
      let updatedDomainList = {
        ...this.state.DomainList,
        [Domain]: deleteDomain,
      };
      this.setState({insertIndex: currentIndex-1});
      this.setState({DomainList: updatedDomainList});
      return currentIndex-1; // fix this
  }

  render() {
    return (
      <div className='DomainSearch'>
        <h3>Construct Module</h3>
        <div className="DomainToolbox">
          <div className="DomainButtonList">
            {this.getAllDomains().map(DomainButton => (
              <Button className='addDomainButton' onClick={ () => {this.insertDomain(DomainButton.domainName)} }>
                {DomainButton.domainName}
              </Button>
              ))
            }
          </div>
          <div className="DomainSandbox">
            {this.getPresentDomains().map(DomainDiv => (
                <Draggable handle='.handle' bounds='parent' axis='x' grid='[36px, 36px]'>
                  <div className={"Domain " + DomainDiv.domainName}>
                    {DomainDiv.domainName}
                    <div className="handle"> -> </div>
                    <div className="deleteIcon" onClick={ () => {this.deleteDomain(DomainDiv.domainName)} }> D </div>
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